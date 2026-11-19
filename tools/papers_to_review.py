#!/usr/bin/env python
"""Render the aggregated paper candidates into a review chunk the owner can sign off.

    py tools\\papers_aggregate.py
    py tools\\papers_to_review.py --chunk review\\09_papers_batches_0_3
    py tools\\apply_review.py review\\09_papers_batches_0_3            (dry run)
    py tools\\apply_review.py review\\09_papers_batches_0_3 --apply    (owner, after deciding)

This writes ONLY into the chunk folder. It never touches data/. Every row arrives with an empty
`decision`, so a chunk that is rendered and then applied without being read changes nothing --
apply_review.py counts a blank decision as undecided and skips it.

What it produces
----------------
entry_proposals.csv   adds readouts to EXISTING entries. One row per entry, not per readout:
                      apply_review.py sets a whole cell, so two rows touching the same
                      (entry_id, readouts) would make the second a conflict against the first.
                      All of an entry's new readouts are therefore merged into one proposed_value.
new_entries.csv       draft rows for functions no playbook entry covers. Ids are LIT1, LIT2 ... so
                      a literature-derived entry is recognisable forever and cannot collide with a
                      hand-curated id.
term_proposals.csv    vocabulary, carrying the scripts/term_table.csv columns so an accepted row
                      can be written there. Applied by apply_review, but accepting a term means
                      the m2 relabel has to be re-run -- see the warning in REVIEW.md.
REVIEW.md             what each file is, how the numbers were reached, and what to check.

The ratings are left blank on purpose, and that blocks the row
--------------------------------------------------------------
screening_suitability and biological_relevance are expert judgements about whether an assay is
worth running, and no number of papers settles them. A paper proves an assay is used, which is a
different claim, so this file will not invent a rating.

Those two fields are also required: tools/build_db.py refuses to build docs/assay_db.json at all if
either is blank, and one bad row takes down the whole site build. apply_review.py therefore refuses
a new entry whose ratings are blank, naming the row. That is deliberate -- a literature-derived
entry cannot reach the site until a person has rated it, which is the correct gate. Accepting a row
in new_entries.csv means filling both fields in the same pass.

Category, role and scope ARE filled, because they have defensible answers from the record itself
and the row cannot be grouped or validated without them. Check the inference rule below.
"""
import argparse
import csv
import os
import re
import sys
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
PAPERS = os.path.join(ROOT, "build_output", "papers")

REVIEW_COLS = ["evidence", "confidence", "decision", "your_note"]
TODAY = datetime.now().strftime("%Y-%m-%d")


def log(msg):
    print(msg)
    sys.stdout.flush()


def read_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))


def header_of(path):
    with open(path, newline="", encoding="utf-8-sig") as h:
        return csv.DictReader(h).fieldnames or []


def split(text):
    return [t.strip() for t in re.split(r"[;|]", text or "") if t.strip()]


def write_rows(path, cols, rows):
    with open(path, "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def num(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return 0


# Fields the owner fills in by hand on a new_entries row and must never lose to a regeneration.
# The ratings are here because apply_review refuses a row without them, so they cost real effort.
OWNER_FIELDS = ["decision", "your_note", "screening_suitability", "biological_relevance",
                "development_priority", "rating_notes", "assay_name", "assay_category", "role",
                "measures", "conditions", "reagents", "readouts"]


def carry_over(rows, previous, key_of, fields):
    """Copy what the owner wrote on a matching earlier row into the freshly generated one.

    The chunk is regenerated every time a batch lands, so without this a regeneration silently
    discards every decision already recorded -- which would make "re-render after each batch" a
    destructive instruction. A generated field is only kept from the old row where the owner
    actually changed it, so new evidence still updates the evidence and count columns.
    """
    old = {key_of(p): p for p in previous}
    kept = 0
    for r in rows:
        p = old.get(key_of(r))
        if not p:
            continue
        touched = False
        for f in fields:
            was = (p.get(f) or "").strip()
            if was and was != (r.get(f) or "").strip():
                r[f] = was
                touched = touched or f in ("decision", "your_note") or bool(was)
        if touched:
            kept += 1
    return kept


# ---------------------------------------------------------------- category and role

def infer_category(function_label, scope):
    """The one inference this file makes. Stated in the docstring and in REVIEW.md so it can be
    checked at a glance rather than trusted."""
    f = function_label.lower()
    if "target engagement" in f or "binding affinity" in f:
        return "biophysical_binding"
    if scope == "cell-based":
        return "cell_functional"
    return "biochemical_hts"


def infer_role(category):
    return "confirmation" if category == "biophysical_binding" else "primary"


def measures_text(function_label, scope):
    return "%s, measured in %s format" % (function_label, scope)


# ---------------------------------------------------------------- proposals

def build_entry_proposals(readout_rows, assays):
    """Readout additions to existing entries, one row per entry.

    A candidate is keyed by (class, scope, readout) and carries `entries_seen`: the entries that the
    underlying records actually matched. Those are the entries entitled to the readout -- a readout
    seen in a kinase-activity record must not land on the selectivity-panel row just because both
    are protein_kinase. A candidate with no entries_seen has nowhere to go and is reported instead,
    because it means the playbook has no row for that context at all.
    """
    by_id = {a["entry_id"]: a for a in assays}
    want, support, orphans = {}, {}, []
    for r in readout_rows:
        targets = [e for e in split(r.get("entries_seen")) if e in by_id]
        if not targets:
            orphans.append(r)
            continue
        for e in targets:
            cur = split(by_id[e].get("readouts"))
            if r["readout"].lower() in {c.lower() for c in cur}:
                continue
            want.setdefault(e, [])
            if r["readout"] not in want[e]:
                want[e].append(r["readout"])
            s = support.setdefault(e, {"papers": 0, "pmcids": set(), "quote": "", "classes": set()})
            s["papers"] = max(s["papers"], num(r.get("n_papers")))
            s["pmcids"].update(split(r.get("pmcids")))
            s["quote"] = s["quote"] or (r.get("quote") or "")
            s["classes"].add(r.get("class_id") or "")

    out = []
    for e, readouts in sorted(want.items()):
        a, s = by_id[e], support[e]
        cur = (a.get("readouts") or "").strip()
        merged = "; ".join(split(cur) + readouts)
        ev = ("%d paper(s) report %s for %s; not listed today. %s. PMC: %s"
              % (s["papers"], ", ".join(readouts), ", ".join(sorted(c for c in s["classes"] if c)),
                 (s["quote"] or "")[:160], ", ".join(sorted(s["pmcids"])[:6])))
        out.append({"entry_id": e, "field": "readouts", "current_value": cur,
                    "proposed_value": merged, "evidence": ev,
                    "confidence": "high" if s["papers"] >= 4 else "medium",
                    "decision": "", "your_note": ""})
    return out, orphans


def entry_key(classes, scope, function_label):
    """What makes a drafted entry the same entry across runs.

    Ids must not be positional. The crawl grows, so the candidate list is re-sorted on every run,
    and if LIT1 were simply "the first row" it would name a different assay next time -- silently
    reassigning whatever decision the owner had already recorded against it. Keying on what the
    entry IS keeps LIT1 the kinase cellular-engagement row for good.
    """
    return ("%s|%s|%s" % (classes, scope, function_label)).strip().lower()


LEDGER = "decided.csv"
LEDGER_COLS = ["key", "entry_id"] + OWNER_FIELDS


def load_ledger(chunk, previous):
    """Everything this chunk must not forget: the key -> LIT id mapping and what the owner wrote.

    Held in its own file rather than inferred from new_entries.csv, because that file only lists the
    candidates of the current run. A candidate that drops below the threshold -- a narrower
    `--batches`, a stricter `--min-papers`, or simply a function whose papers were reclassified --
    would otherwise take both its id and its decision out of circulation with it, and the next
    regeneration would hand that id to a different function. Rows here are never removed, only
    added and updated, so a candidate that comes back comes back decided.
    """
    rows = read_rows(os.path.join(chunk, LEDGER))
    led = {r["key"]: r for r in rows if r.get("key")}
    # new_entries.csv is the file the owner actually edits, so it wins over the ledger's own copy:
    # the ledger was written by the previous run, before those edits existed. Only non-empty values
    # are taken, so clearing a cell in the ledger is not undone by a blank in new_entries.
    for p in previous:
        if not p.get("entry_id"):
            continue
        k = entry_key(p.get("classes"), p.get("scope"), p.get("evidence_labels"))
        rec = led.setdefault(k, {"key": k, "entry_id": p["entry_id"]})
        rec["entry_id"] = rec.get("entry_id") or p["entry_id"]
        for f in OWNER_FIELDS:
            if (p.get(f) or "").strip():
                rec[f] = p[f]
    return led


def save_ledger(chunk, ledger, new_entries):
    """Fold what the current rows carry back into the ledger, then write it."""
    for r in new_entries:
        k = r.get("_key") or ""
        if not k:
            continue
        rec = ledger.setdefault(k, {"key": k, "entry_id": r["entry_id"]})
        rec["entry_id"] = r["entry_id"]
        for f in OWNER_FIELDS:
            if (r.get(f) or "").strip():
                rec[f] = r[f]
    write_rows(os.path.join(chunk, LEDGER), LEDGER_COLS,
               sorted(ledger.values(), key=lambda r: num(re.sub(r"\D", "", r.get("entry_id") or "0"))))


def build_new_entries(entry_rows, assays, assays_header, ledger):
    """Draft entries for functions no row covers.

    Ids and the owner's own edits come from `ledger` for any candidate drafted before, so
    regenerating after a new batch neither renumbers nor un-decides a row. New keys are added to
    `ledger` in place.
    """
    taken = {a["entry_id"] for a in assays}
    taken |= {r.get("entry_id") for r in ledger.values() if r.get("entry_id")}
    n = 1 + max([num(m.group(1)) for m in
                 (re.match(r"^LIT(\d+)$", t) for t in taken) if m] or [0])
    out = []
    for r in entry_rows:
        cat = infer_category(r["function_label"], r["scope"])
        key = entry_key(r["class_id"], r["scope"], r["function_label"])
        rec = ledger.get(key)
        eid = (rec or {}).get("entry_id")
        if not eid:
            eid = "LIT%d" % n
            while eid in taken:
                n += 1
                eid = "LIT%d" % n
            n += 1
        taken.add(eid)
        row = {c: "" for c in assays_header}
        row.update({
            "entry_id": eid,
            "assay_name": "%s (%s)" % (r["function_label"], r["scope"]),
            "classes": r["class_id"],
            "assay_category": cat,
            "scope": r["scope"],
            "role": infer_role(cat),
            "measures": measures_text(r["function_label"], r["scope"]),
            "readouts": r.get("readouts_seen") or "",
            "evidence_function": r["function_label"],
            "evidence_labels": r["function_label"],
            "rationale": ("%s paper(s) across %s protein(s) of this class report this function and no "
                          "playbook entry covers it in %s format."
                          % (r.get("n_papers"), r.get("n_proteins"), r["scope"])),
            "status": "draft from literature %s; needs expert review" % TODAY,
            "reviewer": "",
            "reviewed_on": "",
        })
        row.update({
            "evidence": "Proteins: %s. PMC: %s. %s" % (r.get("proteins"), r.get("pmcids"),
                                                       (r.get("quote") or "")[:160]),
            "confidence": "high" if num(r.get("n_proteins")) >= 3 and num(r.get("n_papers")) >= 3
                          else "medium" if num(r.get("n_proteins")) >= 2 else "low",
            "decision": "", "your_note": "",
        })
        # Anything the owner wrote on this candidate before wins over the generated draft. The
        # evidence and count columns are NOT in OWNER_FIELDS, so they still refresh with new papers.
        for f in OWNER_FIELDS:
            if rec and (rec.get(f) or "").strip():
                row[f] = rec[f]
        row["_key"] = key
        out.append(row)
    return out


def build_system_proposals(system_rows, assays):
    """`system` additions to existing entries, one row per entry, like the readouts.

    `system` is a column data/assays.csv does not have until the first of these is accepted;
    apply_review.py adds it then, because it is listed in that file's NEW_COLUMNS.
    """
    by_id = {a["entry_id"]: a for a in assays}
    want, support = {}, {}
    for r in system_rows:
        for e in [x for x in split(r.get("entries_seen")) if x in by_id]:
            cur = split(by_id[e].get("system"))
            if r["system"].lower() in {c.lower() for c in cur}:
                continue
            want.setdefault(e, [])
            if r["system"] not in want[e]:
                want[e].append(r["system"])
            s = support.setdefault(e, {"papers": 0, "pmcids": set(), "quote": ""})
            s["papers"] = max(s["papers"], num(r.get("n_papers")))
            s["pmcids"].update(split(r.get("pmcids")))
            s["quote"] = s["quote"] or (r.get("quote") or "")
    out = []
    for e, systems in sorted(want.items()):
        cur = (by_id[e].get("system") or "").strip()
        s = support[e]
        out.append({"entry_id": e, "field": "system", "current_value": cur,
                    "proposed_value": "; ".join(split(cur) + systems),
                    "evidence": "%d paper(s) name %s as the preparation. %s. PMC: %s"
                                % (s["papers"], ", ".join(systems), (s["quote"] or "")[:160],
                                   ", ".join(sorted(s["pmcids"])[:6])),
                    "confidence": "high" if s["papers"] >= 4 else "medium",
                    "decision": "", "your_note": ""})
    return out


# scripts/term_table.csv, the labeller's vocabulary. A term proposal carries these columns so an
# accepted row can be written straight into that file, the way new_entries.csv carries the
# assays.csv columns. Only `technology` and `function` exist there as axes: a `system` or
# `counter_screen` proposal is rendered with a blank axis, and apply_review refuses it by name
# rather than guessing where it belongs.
TERM_TABLE_COLS = ["axis", "term", "synonyms", "examples", "suppresses", "applies_to", "notes"]
TERM_AXES = {"technology", "function"}


def build_term_proposals(term_rows):
    out = []
    for t in term_rows:
        kind = (t.get("kind") or "").strip().lower()
        variants = t.get("variants") or ""
        first = next((v.strip() for v in variants.split("|") if v.strip()), "")
        # The labeller reads ChEMBL descriptions, which phrase things differently from a methods
        # section. label_synonyms holds that wording and comes first, because it is the form the
        # term will actually be matched by; the paper verbatims follow as extra cover.
        syn, seen_syn = [], set()
        for s in ((t.get("label_synonyms") or "") + "|" + variants).split("|"):
            s = s.strip()
            if s and s.lower() not in seen_syn:
                seen_syn.add(s.lower())
                syn.append(s)
        out.append({
            # the row as it would be written into scripts/term_table.csv
            "axis": kind if kind in TERM_AXES else "",
            "term": t.get("proposed_term"),
            "synonyms": "|".join(syn),
            "examples": first or (t.get("definition") or "")[:120],
            "suppresses": t.get("label_suppresses") or "",
            "applies_to": "any",
            "notes": "from the paper crawl, %s: %s" % (t.get("term_id") or "?",
                                                       (t.get("definition") or "")[:200]),
            # the evidence behind it, for the decision
            "term_id": t.get("term_id"), "proposed_term": t.get("proposed_term"), "kind": kind,
            "nearest_controlled_value": t.get("nearest_controlled_value"),
            "definition": t.get("definition"), "variants": variants,
            "n_papers": t.get("n_papers"), "n_proteins": t.get("n_proteins"),
            "proteins": t.get("proteins"), "pmcids": t.get("pmcids"),
            "evidence": (t.get("quote") or "")[:200],
            "confidence": "high" if num(t.get("n_papers")) >= 4 else "medium",
            "decision": "", "your_note": ""})
    return out


# ---------------------------------------------------------------- the readme

def review_md(chunk, props, orphans, new_entries, terms, min_papers, problems, systems_waiting):
    name = os.path.basename(chunk.rstrip("\\/"))
    L = ["# %s" % name, "",
         "Paper-derived proposals, generated %s by `tools/papers_to_review.py` from"
         " `build_output/papers/*_candidates.csv` (written by `tools/papers_aggregate.py`)." % TODAY,
         "",
         "Every row has an empty `decision`. Set it to `accept`, `reject`, or `edit` (and put the"
         " value you want in `your_note`), then:", "",
         "```", "py tools\\apply_review.py %s            # dry run, shows what would change" % name,
         "py tools\\apply_review.py %s --apply" % name, "```", "",
         "## How the numbers were reached", "",
         "- A candidate needs **%d distinct papers**. One paper is an anecdote." % min_papers,
         "- Candidates are aggregated **per function class, not per protein**. `data/assays.csv` is"
         " keyed by class, so a kinase paper improves `protein_kinase` and every kinase inherits it."
         " `n_proteins` tells you how many distinct proteins stand behind a row: 1 protein at 3"
         " papers is weaker than 3 proteins at 3 papers, and the CSV carries both.",
         "- Records with `confidence: low` were excluded from counting.",
         "- **Check any row where `n_proteins` is 1.** Two papers about one protein is the weakest"
         " thing above the threshold, and it is also where a misclassified protein shows up: the"
         " class a candidate is filed under comes from the survey's classification, which is"
         " sometimes wrong. `solute_carrier` backed only by AKT1 is a classification error, not a"
         " transporter finding. Such rows are marked `confidence: low`.",
         "- Nothing here enters the campaign counts. Paper evidence is a separate tier"
         " (`data/literature_stats.csv`, written by `scripts/literature_stats.py`).", "",
         "## Files", ""]

    L += ["### entry_proposals.csv - %d row(s)" % len(props), "",
          "Adds readouts to entries that already exist. **One row per entry, not per readout**:"
          " `apply_review.py` writes a whole cell, so the `proposed_value` is the entry's current"
          " `readouts` plus every new one. Rejecting the row rejects all of them; to take some and"
          " not others, set `decision = edit` and write the list you want in `your_note`.", ""]
    if props:
        L += ["| entry | adds | papers |", "|---|---|---|"]
        for p in props[:30]:
            added = [r for r in split(p["proposed_value"]) if r not in split(p["current_value"])]
            L.append("| %s | %s | %s |" % (p["entry_id"], ", ".join(added),
                                           re.search(r"^(\d+)", p["evidence"]).group(1)
                                           if re.search(r"^(\d+)", p["evidence"]) else "?"))
        L.append("")

    L += ["### new_entries.csv - %d row(s)" % len(new_entries), "",
          "Functions that **no** playbook entry covers for that class in that format. These are the"
          " expensive gaps: not a missing readout on an existing row, a missing row.", "",
          "**To accept a row you must fill `screening_suitability` and `biological_relevance`**"
          " (one of very_low, low, moderate, high, very_high). They are blank because papers prove"
          " an assay is *used*, not that it is worth running, and this generator will not invent a"
          " rating. They are also required: `build_db.py` refuses to build the site if either is"
          " blank, so `apply_review.py` refuses the row first and names it. A literature-derived"
          " entry cannot reach the site until a person has rated it.", "",
          "`assay_category` and `role` **are** inferred and are the two fields to check: a function"
          " naming target engagement or binding affinity becomes `biophysical_binding` /"
          " `confirmation`, anything else cell-based becomes `cell_functional` / `primary`, and the"
          " rest `biochemical_hts` / `primary`.", ""]
    if new_entries:
        L += ["| id | class | scope | function | papers | proteins |", "|---|---|---|---|---|---|"]
        for r in new_entries[:40]:
            prot = re.search(r"Proteins: ([^.]*)\.", r["evidence"])
            L.append("| %s | %s | %s | %s | %s | %s |" % (
                r["entry_id"], r["classes"], r["scope"], r["evidence_labels"],
                (re.search(r"^(\d+) paper", r["rationale"]).group(1)
                 if re.search(r"^(\d+) paper", r["rationale"]) else "?"),
                (prot.group(1)[:40] if prot else "")))
        L.append("")

    L += ["### term_proposals.csv - %d row(s)" % len(terms), "",
          "`apply_review.py` writes an accepted row into `scripts/term_table.csv`. **Accepting a"
          " term is still the expensive decision in this repository:** a new technology term"
          " changes how every campaign is labelled, so the full m2 relabel has to be re-run before"
          " the usage statistics mean anything. Batch every term you intend to accept into one"
          " pass; do not take them a few at a time.", "",
          "A `system` or `counter_screen` term has no axis in `term_table.csv` and is refused by"
          " name; decide where those belong before accepting them.", ""]
    if terms:
        L += ["| term | kind | nearest existing | papers |", "|---|---|---|---|"]
        for t in terms[:40]:
            L.append("| %s | %s | %s | %s |" % (t["proposed_term"], t["kind"],
                                                t["nearest_controlled_value"] or "-", t["n_papers"]))
        L.append("")

    if orphans:
        L += ["## Readouts with nowhere to go - %d" % len(orphans), "",
              "These readouts are above the threshold but every record carrying them matched **no"
              " entry**, so there is no cell to add them to. They are not dropped: the same records"
              " produced the `new_entries.csv` rows above, and once those entries exist these"
              " readouts belong on them. Re-run the aggregator after accepting new entries.", "",
              "| class | scope | readout | papers |", "|---|---|---|---|"]
        for r in orphans[:30]:
            L.append("| %s | %s | %s | %s |" % (r["class_id"], r["scope"], r["readout"],
                                                r["n_papers"]))
        L.append("")

    if systems_waiting:
        L += ["## Systems and counter-screens, not yet proposable - %d" % len(systems_waiting), "",
              "The `system` column (the biological preparation an assay runs in: a Ba/F3 line,"
              " HepaRG, a donor microsome panel, an organoid) now exists end to end -- a `system`"
              " field proposal creates it in `data/assays.csv`, `build_db.py` publishes it and the"
              " site shows a **Systems used** line. Nothing can fill it yet.", "",
              "**Every value below rests on exactly one paper.** That is not a shortage of papers,"
              " it is the shape of the evidence: a methods section names its preparation once, and"
              " two papers on the same target rarely choose the same one. Corroboration therefore"
              " needs more papers *in the same class*, not more papers overall -- which is the"
              " strongest argument for picking the next batches by empty class.", "",
              "| class | scope | value | kind | protein | paper |", "|---|---|---|---|---|---|"]
        for r in systems_waiting[:40]:
            L.append("| %s | %s | %s | %s | %s | %s |" % (
                r.get("class_id"), r.get("scope"), r.get("system"), r.get("kind"),
                r.get("proteins"), r.get("pmcids")))
        L.append("")

    if problems:
        L += ["## Records that contradict themselves - %d" % len(problems), "",
              "Nothing is proposed from these. A function label that says *cellular* cannot be a"
              " cell-free gap, so the record is wrong rather than the playbook. All of them are"
              " lysate methods -- DARTS, CETSA on lysate -- recorded as `cell-free` while keeping"
              " the cellular-engagement function. Fix them in the extraction if they matter;"
              " proposing an entry for them would write the contradiction into the database.", "",
              "| function | protein | paper | readout |", "|---|---|---|---|"]
        for r in problems[:20]:
            L.append("| %s | %s | %s | %s |" % (r.get("function_label"), r.get("gene"),
                                                r.get("pmcid"), r.get("technology_labels")))
        L.append("")
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chunk", required=True, help="chunk folder to write, e.g. review\\09_papers")
    ap.add_argument("--min-papers", type=int, default=2, help="echoed into REVIEW.md")
    ap.add_argument("--dry-run", action="store_true", help="print the counts, write nothing")
    args = ap.parse_args(argv)

    assays_path = os.path.join(DATA, "assays.csv")
    assays, assays_header = read_rows(assays_path), header_of(assays_path)
    if not assays:
        raise SystemExit("data/assays.csv not readable")

    readout_rows = read_rows(os.path.join(PAPERS, "readout_candidates.csv"))
    entry_rows = read_rows(os.path.join(PAPERS, "new_entry_candidates.csv"))
    term_rows = read_rows(os.path.join(PAPERS, "term_candidates.csv"))
    problems = read_rows(os.path.join(PAPERS, "record_problems.csv"))
    system_rows = read_rows(os.path.join(PAPERS, "system_candidates.csv"))
    systems_waiting = read_rows(os.path.join(PAPERS, "system_candidates_below_threshold.csv"))
    if not (readout_rows or entry_rows or term_rows):
        raise SystemExit("no candidate files: run py tools\\papers_aggregate.py first")

    # A chunk is regenerated as batches land, so read what is already there before overwriting it.
    prev_props = read_rows(os.path.join(args.chunk, "entry_proposals.csv"))
    prev_entries = read_rows(os.path.join(args.chunk, "new_entries.csv"))
    prev_terms = read_rows(os.path.join(args.chunk, "term_proposals.csv"))

    ledger = load_ledger(args.chunk, prev_entries)
    props, orphans = build_entry_proposals(readout_rows, assays)
    props += build_system_proposals(system_rows, assays)
    new_entries = build_new_entries(entry_rows, assays, assays_header, ledger)
    terms = build_term_proposals(term_rows)

    # new_entries already took its decisions from the ledger, which outlives a row dropping out.
    kept = sum(1 for r in new_entries if (r.get("decision") or "").strip())
    kept += carry_over(props, prev_props,
                       lambda r: (r.get("entry_id"), r.get("field")), ["decision", "your_note"])
    kept += carry_over(terms, prev_terms, lambda r: r.get("term_id"), ["decision", "your_note"])

    log("From %d readout, %d entry and %d term candidate(s):"
        % (len(readout_rows), len(entry_rows), len(term_rows)))
    log("  entry_proposals.csv   %4d row(s)  (readouts onto existing entries)" % len(props))
    log("  new_entries.csv       %4d row(s)  (no playbook row covers this)" % len(new_entries))
    log("  term_proposals.csv    %4d row(s)  (-> scripts/term_table.csv; forces an m2 relabel)" % len(terms))
    if orphans:
        log("  %d readout candidate(s) matched no entry; listed in REVIEW.md, not dropped."
            % len(orphans))
    if kept:
        log("  %d row(s) kept what you had already written (decisions, ratings, wording)." % kept)
    if args.dry_run:
        log("\nDry run: nothing written.")
        return 0

    os.makedirs(args.chunk, exist_ok=True)
    save_ledger(args.chunk, ledger, new_entries)
    write_rows(os.path.join(args.chunk, "entry_proposals.csv"),
               ["entry_id", "field", "current_value", "proposed_value"] + REVIEW_COLS, props)
    write_rows(os.path.join(args.chunk, "new_entries.csv"),
               assays_header + REVIEW_COLS, new_entries)
    write_rows(os.path.join(args.chunk, "term_proposals.csv"),
               TERM_TABLE_COLS +
               ["term_id", "proposed_term", "kind", "nearest_controlled_value", "definition",
                "variants", "n_papers", "n_proteins", "proteins", "pmcids"] + REVIEW_COLS, terms)
    with open(os.path.join(args.chunk, "REVIEW.md"), "w", encoding="utf-8") as h:
        h.write(review_md(args.chunk, props, orphans, new_entries, terms, args.min_papers,
                          problems, systems_waiting))

    log("\nWritten to %s (%s)." % (args.chunk, "decisions you had already made were kept"
                                   if kept else "every decision blank"))
    if new_entries:
        log("Note: accepting a new_entries.csv row needs screening_suitability and")
        log("      biological_relevance filled in too; apply_review refuses the row otherwise.")
    log("Next: py tools\\apply_review.py %s" % args.chunk)
    return 0


if __name__ == "__main__":
    sys.exit(main())
