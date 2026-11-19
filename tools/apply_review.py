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
  term_proposals.csv   NEW vocabulary terms (term_table.csv columns + review columns)    -> scripts/term_table.csv
                       Accepting one changes how every campaign is labelled: the m2 relabel has to
                       be re-run afterwards, so batch the terms you mean to accept into one pass.

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
          "classes": ("function_classes.csv", "class_id"),
          "literature": ("literature_stats.csv", "lit_id"),
          "terms": ("term_table.csv", "term")}
# term_table.csv is the labeller's vocabulary and lives next to the script that reads it, not in
# data/. It is the one table here whose key is the human-readable value itself, so a duplicate
# check is a check on the term's wording.
OTHER_DIRS = {"terms": "scripts"}
# A term may only be filed under an axis that m2_label.py understands.
TERM_AXES = {"technology", "function"}
# Tables that may not exist yet. literature_stats.csv is created by the first accepted literature
# chunk; every other table is part of the repository and its absence is a real problem.
CREATABLE = {"literature": ["lit_id", "entry_id", "scope_type", "scope_id", "effect", "readout",
                            "n_papers_total", "n_papers_using", "pmcids", "evidence", "confidence",
                            "source", "retrieved_on", "method_version", "notes"]}
# Columns a chunk may introduce: when a proposal names one that the table does not have yet, the column is added
# (empty for every other row). Any other unknown column is still an error.
NEW_COLUMNS = {"assays": ["qc_metric", "dmso_pct", "interference", "system"]}
# Mirrors tools/build_db.py, which is the source of truth. A new assay row that breaks these is
# refused here rather than after it has been written: build_db then REFUSES TO BUILD THE SITE AT
# ALL, so one bad row in a chunk takes down the whole output, and the error surfaces a step later
# than the decision that caused it. Blank is not allowed for either level.
ASSAY_LEVELS = {"very_low", "low", "moderate", "high", "very_high"}
ASSAY_ROLES = {"primary", "confirmation", "counterscreen"}
ASSAY_CATEGORIES = {"biochemical_hts", "cell_functional", "biophysical_binding",
                    "orthogonal_confirmation", "mechanistic_counterscreen"}
ASSAY_SCOPES = {"cell-free", "cell-based"}
# (file in the chunk, kind, table). Classes first, then entries, then kits, so new rows can be referred to.
OPERATIONS = [("new_classes.csv", "new", "classes"), ("class_updates.csv", "field", "classes"),
              ("new_entries.csv", "new", "assays"), ("entry_proposals.csv", "field", "assays"),
              ("kit_proposals.csv", "new", "kits"), ("kit_updates.csv", "field", "kits"),
              ("literature_proposals.csv", "new", "literature"),
              ("term_proposals.csv", "new", "terms")]


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
    if proposals and key not in proposals[0]:
        # A chunk rendered by an older papers_to_review.py: the file predates the columns of the
        # table it now feeds. Say so once, rather than raising KeyError on the first row.
        return [("error", "-", "this file has no '%s' column, so it was written for an older "
                               "version of this tool: re-render the chunk with "
                               "tools\\papers_to_review.py" % key)]
    for p in proposals:
        rid, d = (p.get(key) or "").strip(), decision(p)
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
        # most tables live in data/; term_table.csv lives beside the script that reads it. A test
        # pointing --data at its own copy still wins, so the fixtures stay self-contained.
        path = os.path.join(args.data, fname)
        if name in OTHER_DIRS and not os.path.exists(path):
            path = os.path.join(ROOT, OTHER_DIRS[name], fname)
        if os.path.exists(path):
            head, rows, lay = read_table(path)
        elif name in CREATABLE and os.path.exists(os.path.join(args.chunk,
                                                               "%s_proposals.csv" % name)):
            # First accepted chunk of this kind: start the file, matching the house layout.
            head, rows, lay = CREATABLE[name], [], {"bom": True, "eol": "\r\n"}
            print("Creating data/%s" % fname)
        else:
            continue
        tables[name] = {"path": path, "head": head, "rows": rows, "lay": lay, "key": key,
                        "changed": False, "new_file": not os.path.exists(path)}

    # ids that exist after this chunk (for checking references to rows added in the same chunk)
    def ids(table):
        return {r[TABLES[table][1]] for r in tables[table]["rows"]} if table in tables else set()
    class_ids, entry_ids = ids("classes"), ids("assays")

    def check_entry(row):
        """Refuse a new assay row that data/assays.csv would accept but build_db.py would not.

        Each rule applies only if the table actually has that column: assays.csv has them all, but
        a caller pointing --data at a reduced table should not be told a column is wrong when there
        is no such column to be wrong about."""
        unknown = [c for c in ids_in(row.get("classes")) if c not in class_ids]
        if unknown and "classes" in tables:
            return "uses unknown classes: " + ", ".join(unknown)
        head = set(tables.get("assays", {}).get("head") or [])
        blank = [c for c in ("screening_suitability", "biological_relevance")
                 if c in head and (row.get(c) or "") not in ASSAY_LEVELS]
        if blank:
            return ("%s must be one of %s; blank is not allowed because build_db.py would then "
                    "refuse to build the site" % (", ".join(blank), ", ".join(sorted(ASSAY_LEVELS))))
        # "primary or confirmation" is in use, so only the first word is checked, as build_db does.
        if "role" in head and (row.get("role") or "").replace(" or ", ";").split(";")[0].strip().lower() \
                not in ASSAY_ROLES:
            return "role must be one of %s, got '%s'" % (", ".join(sorted(ASSAY_ROLES)), row.get("role"))
        if "assay_category" in head and (row.get("assay_category") or "") not in ASSAY_CATEGORIES:
            return "assay_category must be one of %s, got '%s'" % (
                ", ".join(sorted(ASSAY_CATEGORIES)), row.get("assay_category"))
        if "scope" in head and (row.get("scope") or "") not in ASSAY_SCOPES:
            return "scope must be one of %s, got '%s'" % (", ".join(sorted(ASSAY_SCOPES)),
                                                          row.get("scope"))
        return ""

    def check_kit(row):
        unknown = [e for e in ids_in(row.get("applies_to_entries")) if e not in entry_ids]
        return "applies to unknown entries: " + ", ".join(unknown) if unknown else ""

    def check_term(row):
        """A vocabulary row that m2_label.py has to be able to use.

        The axis is the whole of it: m2_label builds one matcher per axis, so a row filed under
        anything else is silently never matched. `system` and `counter_screen` proposals arrive
        with a blank axis on purpose -- where they belong is an open question, and guessing here
        would bury it."""
        axis = (row.get("axis") or "").strip().lower()
        if not axis:
            return ("no axis: a '%s' term has no home in term_table.csv yet"
                    % ((row.get("kind") or "?").strip() or "?"))
        if axis not in TERM_AXES:
            return "axis must be one of %s, got '%s'" % (", ".join(sorted(TERM_AXES)), axis)
        if not (row.get("term") or "").strip():
            return "no term"
        if (row.get("applies_to") or "any").strip().lower() not in ("any", "cell-based"):
            return "applies_to must be 'any' or 'cell-based', got '%s'" % row.get("applies_to")
        return ""

    def check_literature(row):
        """A literature row keys onto a usage row, so its entry must exist and its counts must be
        numbers in the right order; build_database would otherwise silently skip it later."""
        if entry_ids and row.get("entry_id") not in entry_ids:
            return "unknown entry '%s'" % row.get("entry_id")
        if row.get("scope_type") not in ("gene", "function_class"):
            return "scope_type must be gene or function_class"
        # "*" means every effect view, for a finding with no direction (see literature_stats.py).
        if (row.get("effect") or "") not in ("*", "inhibition", "activation",
                                             "inhibition_labeled", "activation_with_unspecified"):
            return "effect must be an effect view or '*' for all of them"
        try:
            total, using = int(row["n_papers_total"]), int(row["n_papers_using"])
        except (KeyError, ValueError):
            return "n_papers_total and n_papers_using must be whole numbers"
        if using > total:
            return "n_papers_using (%d) is larger than n_papers_total (%d)" % (using, total)
        return ""

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
            check = {"assays": check_entry, "kits": check_kit,
                     "literature": check_literature, "terms": check_term}.get(table)
            plan = plan_new_rows(props, t["head"], t["rows"], t["key"], check)
            new_ids = {rid for st, rid, _ in plan if st in ("change", "already")}
            if table == "classes":
                class_ids |= new_ids
            if table == "assays":
                entry_ids |= new_ids
        plans.append((fname, kind, table, plan))

    for fname, kind, table, plan in plans:
        print("%s -> %s/%s:" % (fname, OTHER_DIRS.get(table, "data"), TABLES[table][0]))
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
                if st == "change":
                    label = (info.get("product") or info.get("assay_name") or info.get("label")
                             or ("%s %s %s/%s papers" % (info.get("scope_id", ""),
                                                         info.get("readout") or "(whole entry)",
                                                         info.get("n_papers_using", ""),
                                                         info.get("n_papers_total", ""))
                                 if info.get("scope_id") else ""))
                else:
                    label = info
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
            if not t.get("new_file"):       # nothing to back up the first time a table is created
                shutil.copy2(t["path"], os.path.join(bdir, "%s_%s.csv" % (base, stamp)))
            write_table(t["path"], t["head"], t["rows"], t["lay"])
    print("Written. Backups: data\\backup\\*_%s.csv" % stamp)
    print("Next: py scripts\\fetch_citations.py   then   py tools\\build_db.py")
    if tables.get("terms", {}).get("changed"):
        print("scripts\\term_table.csv changed: re-run the m2 relabel before trusting the usage")
        print("statistics, because every campaign's labels were produced by the old vocabulary.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
