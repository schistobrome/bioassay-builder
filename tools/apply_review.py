#!/usr/bin/env python3
"""Copy the changes you ACCEPTED in a review chunk into the curated files in data/.

    py tools\\apply_review.py review\\01_ion_channel            (dry run: only shows what would change)
    py tools\\apply_review.py review\\01_ion_channel --apply    (backs up the files it changes, then writes them)

What it reads from the chunk folder (every file is optional; see review/README.md):
  entry_proposals.csv  one field of one EXISTING assay entry per row -> data/assays.csv
  new_entries.csv      NEW assay entries (assays.csv columns + review columns)          -> data/assays.csv
  kit_proposals.csv    NEW kits (kits.csv columns + review columns)                      -> data/kits.csv
  kit_updates.csv      one field of one EXISTING kit per row                             -> data/kits.csv
  class_updates.csv    one field of one EXISTING function class per row                  -> data/function_classes.csv
  new_classes.csv      NEW function classes (function_classes.csv columns + review columns) -> data/function_classes.csv

Decisions: accept -> use proposed_value (or the row as written); edit -> for field changes use your_note as the new
value, for new rows use the row as you edited it; reject or blank -> skipped.

Safety checks, in plain words:
  - A field is changed only if the file still holds the value the proposal was written against (current_value).
    If you edited that cell since, the row is reported as a CONFLICT and left alone.
  - A new row is added only if its id is not taken. A new kit must apply to entries that exist; a new entry must use
    classes that exist (new ones from the same chunk count).
  - Running it twice does nothing the second time ("already applied").
  - Files keep their exact layout (byte-order mark, Windows line endings), so a diff shows only real changes.
"""
import argparse
import csv
import io
import os
import shutil
import sys
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TABLES = {"assays": ("assays.csv", "entry_id"), "kits": ("kits.csv", "kit_id"),
          "classes": ("function_classes.csv", "class_id")}
# Columns a chunk may introduce: when a proposal names one that the table does not have yet, the column is added
# (empty for every other row). Any other unknown column is still an error.
NEW_COLUMNS = {"assays": ["qc_metric", "dmso_pct", "interference"]}
# (file in the chunk, kind, table). Classes first, then entries, then kits, so new rows can be referred to.
OPERATIONS = [("new_classes.csv", "new", "classes"), ("class_updates.csv", "field", "classes"),
              ("new_entries.csv", "new", "assays"), ("entry_proposals.csv", "field", "assays"),
              ("kit_proposals.csv", "new", "kits"), ("kit_updates.csv", "field", "kits")]


def read_table(path):
    """-> (header, rows as dicts, layout) where layout remembers the BOM and line ending for writing back."""
    with open(path, "rb") as f:
        raw = f.read()
    layout = {"bom": raw.startswith(b"\xef\xbb\xbf"), "eol": "\r\n" if b"\r\n" in raw else "\n"}
    text = raw.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    return reader.fieldnames, list(reader), layout


def write_table(path, header, rows, layout):
    out = io.StringIO()
    w = csv.DictWriter(out, fieldnames=header, lineterminator=layout["eol"])
    w.writeheader()
    w.writerows(rows)
    data = out.getvalue().encode("utf-8")
    with open(path, "wb") as f:
        f.write((b"\xef\xbb\xbf" if layout["bom"] else b"") + data)


def decision(row):
    return (row.get("decision") or "").strip().lower()


def ids_in(cell):
    return [x.strip() for x in (cell or "").replace("|", ";").split(";") if x.strip()]


def plan_entries(proposals, header, rows, key="entry_id"):
    """Field changes to existing rows. -> list of (status, row id, field, old, new).
    status: change / already / conflict / error / skip / undecided"""
    by_id = {a[key]: a for a in rows}
    out = []
    for p in proposals:
        e, field, d = p[key].strip(), p["field"].strip(), decision(p)
        if d in ("", "undecided"):
            out.append(("undecided", e, field, "", "")); continue
        if d == "reject":
            out.append(("skip", e, field, "", "")); continue
        if d not in ("accept", "edit"):
            out.append(("error", e, field, "", "decision must be accept, edit, reject or blank, got '%s'" % d)); continue
        new = (p["your_note"] if d == "edit" else p["proposed_value"]).strip()
        if e not in by_id or field not in header or field == key:
            out.append(("error", e, field, "", "unknown row or column")); continue
        if d == "edit" and not new:
            out.append(("error", e, field, "", "decision is 'edit' but your_note is empty")); continue
        cur = (by_id[e].get(field) or "").strip()
        if cur == new:
            out.append(("already", e, field, cur, new))
        elif cur != (p.get("current_value") or "").strip():
            out.append(("conflict", e, field, cur, new))
        else:
            out.append(("change", e, field, cur, new))
    return out


def plan_new_rows(proposals, header, rows, key, check=None):
    """New rows. -> list of (status, id, row or message). check(row) returns an error message or ''."""
    have = {k[key]: k for k in rows}
    out, seen = [], set()
    for p in proposals:
        rid, d = p[key].strip(), decision(p)
        if d in ("", "undecided"):
            out.append(("undecided", rid, "")); continue
        if d == "reject":
            out.append(("skip", rid, "")); continue
        if d not in ("accept", "edit"):
            out.append(("error", rid, "decision must be accept, edit, reject or blank")); continue
        row = {c: (p.get(c) or "").strip() for c in header}
        problem = check(row) if check else ""
        if problem:
            out.append(("error", rid, problem))
        elif rid in have:
            same = all((have[rid].get(c) or "").strip() == row[c] for c in header)
            out.append(("already" if same else "conflict", rid, "" if same else "%s already used by a different row" % key))
        elif rid in seen:
            out.append(("error", rid, "%s appears twice in the proposals" % key))
        else:
            seen.add(rid)
            out.append(("change", rid, row))
    return out


def plan_kits(proposals, kit_header, kits, entry_ids):
    """New kits; every entry they apply to must exist. (Kept for compatibility.)"""
    def check(row):
        unknown = [e for e in ids_in(row.get("applies_to_entries")) if e not in entry_ids]
        return "applies to unknown entries: " + ", ".join(unknown) if unknown else ""
    return plan_new_rows(proposals, kit_header, kits, "kit_id", check)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("chunk", help="review chunk folder, for example review\\01_ion_channel")
    ap.add_argument("--apply", action="store_true", help="write the changes (default: dry run)")
    ap.add_argument("--data", default=os.path.join(ROOT, "data"), help=argparse.SUPPRESS)  # tests point this elsewhere
    args = ap.parse_args(argv)

    tables = {}
    for name, (fname, key) in TABLES.items():
        path = os.path.join(args.data, fname)
        if os.path.exists(path):
            head, rows, lay = read_table(path)
            tables[name] = {"path": path, "head": head, "rows": rows, "lay": lay, "key": key, "changed": False}

    # ids that exist after this chunk (for checking references to rows added in the same chunk)
    def ids(table):
        return {r[TABLES[table][1]] for r in tables[table]["rows"]} if table in tables else set()
    class_ids, entry_ids = ids("classes"), ids("assays")

    def check_entry(row):
        unknown = [c for c in ids_in(row.get("classes")) if c not in class_ids]
        return "uses unknown classes: " + ", ".join(unknown) if unknown and "classes" in tables else ""

    def check_kit(row):
        unknown = [e for e in ids_in(row.get("applies_to_entries")) if e not in entry_ids]
        return "applies to unknown entries: " + ", ".join(unknown) if unknown else ""

    plans = []          # (label, file, kind, table, plan)
    for fname, kind, table in OPERATIONS:
        path = os.path.join(args.chunk, fname)
        if not os.path.exists(path) or table not in tables:
            continue
        t, props = tables[table], read_table(path)[1]
        for col in NEW_COLUMNS.get(table, []):
            if col not in t["head"] and any((p.get("field") or "").strip() == col or col in p for p in props):
                t["head"] = t["head"] + [col]
                print("New column for data/%s: %s" % (TABLES[table][0], col))
        if kind == "field":
            plan = plan_entries(props, t["head"], t["rows"], key=t["key"])
            if table == "kits":                          # an updated applies_to_entries must still name real entries
                for i, (st, kid, f, old, new) in enumerate(plan):
                    bad = [e for e in ids_in(new) if e not in entry_ids]
                    if st == "change" and f == "applies_to_entries" and bad:
                        plan[i] = ("error", kid, f, old, "applies to unknown entries: " + ", ".join(bad))
        else:
            check = {"assays": check_entry, "kits": check_kit}.get(table)
            plan = plan_new_rows(props, t["head"], t["rows"], t["key"], check)
            new_ids = {rid for st, rid, _ in plan if st in ("change", "already")}
            if table == "classes":
                class_ids |= new_ids
            if table == "assays":
                entry_ids |= new_ids
        plans.append((fname, kind, table, plan))

    for fname, kind, table, plan in plans:
        print("%s -> data/%s:" % (fname, TABLES[table][0]))
        for item in plan:
            if kind == "field":
                st, rid, f, old, new = item
                if st == "change":
                    print("  CHANGE   %-7s %-18s '%s' -> '%s'" % (rid, f, old[:40], new[:80]))
                elif st in ("conflict", "error"):
                    print("  %-8s %-7s %-18s %s" % (st.upper(), rid, f, new if st == "error" else "file now holds '%s'" % old[:60]))
                else:
                    print("  %-8s %-7s %s" % (st, rid, f))
            else:
                st, rid, info = item
                label = (info.get("product") or info.get("assay_name") or info.get("label") or "") if st == "change" else info
                print("  %-8s %-7s %s" % (st.upper() if st in ("change", "conflict", "error") else st, rid, label))

    every = [item for *_, plan in plans for item in plan]
    n_ch = sum(1 for x in every if x[0] == "change")
    n_bad = sum(1 for x in every if x[0] in ("conflict", "error"))
    n_und = sum(1 for x in every if x[0] == "undecided")
    print("\n%d change(s) to make, %d conflict/error(s), %d row(s) without a decision." % (n_ch, n_bad, n_und))
    if not args.apply:
        print("Dry run: nothing written. Add --apply to write." if n_ch else "Nothing to write.")
        return 0
    if not n_ch:
        print("Nothing to write.")
        return 0

    for fname, kind, table, plan in plans:
        t = tables[table]
        if kind == "field":
            by_id = {r[t["key"]]: r for r in t["rows"]}
            for st, rid, f, old, new in plan:
                if st == "change":
                    by_id[rid][f] = new; t["changed"] = True
        else:
            for st, rid, row in plan:
                if st == "change":
                    t["rows"].append(row); t["changed"] = True
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    bdir = os.path.join(args.data, "backup")
    os.makedirs(bdir, exist_ok=True)
    for name, t in tables.items():
        if t["changed"]:
            base = os.path.splitext(TABLES[name][0])[0]
            shutil.copy2(t["path"], os.path.join(bdir, "%s_%s.csv" % (base, stamp)))
            write_table(t["path"], t["head"], t["rows"], t["lay"])
    print("Written. Backups: data\\backup\\*_%s.csv" % stamp)
    print("Next: py scripts\\fetch_citations.py   then   py tools\\build_db.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
