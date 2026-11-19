#!/usr/bin/env python
"""Turn extracted papers into a review chunk the owner can accept, reject or edit.

Where paper evidence goes, and why it goes there
------------------------------------------------
data/usage_stats.csv already has two columns that nothing has ever filled: n_literature_total and
n_literature_using. tools/build_db.py already scores them:

    sc = n_campaigns_using / n_campaigns_classifiable      campaign share
    sl = n_literature_using / n_literature_total           literature share
    score = weighted mean, weights from data/config.json (campaigns 0.7, literature 0.3)

and docs/index.html already renders "...; 60% of 25 papers" on the Function occurrence line. So the
database, the scorer and the site were all built for this. Nothing new has to be invented; the
columns have to be filled.

What must NOT happen is writing those numbers straight into data/usage_stats.csv. That file is
regenerated: build_database.py --merge replaces rows whose method_version starts with "bd-" and
keeps hand-made rows untouched. A hand-made literature row would therefore freeze its campaign
counts at whatever they were the day it was written, and the protein would stop tracking ChEMBL and
PubChem. So paper evidence is an INPUT to generation, held in data/literature_stats.csv, and
build_database.py folds it into the rows it generates. Campaign counts stay live; paper counts ride
along; one rebuild refreshes both.

The flow
--------
    scripts/fetch_papers.py          search, score, gate on subject, split into assay blocks
    (agent or local model)           fill schema/paper_extraction.schema.json per paper
    scripts/literature_stats.py      THIS FILE: paper JSONs -> review/NN_literature/
    owner reviews                    sets decision = accept / reject / edit on each row
    tools/apply_review.py --apply    writes data/literature_stats.csv
    build_database.py --merge        folds it into usage_stats.csv and rebuilds
    tools/build_db.py                docs/assay_db.json, and the site shows it

Every row carries its PMCIDs and a verbatim quote, so a reviewer can check the claim without
opening the pipeline.

Usage
-----
  py scripts\\literature_stats.py --chunk review\\08_literature
  py scripts\\literature_stats.py --chunk review\\08_literature --genes GPBAR1,TAAR1 --dry-run
"""

import argparse
import collections
import csv
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
PAPERS = os.path.join(ROOT, "build_output", "papers")

sys.path.insert(0, HERE)
import target_subject as tsub       # noqa: E402

# The columns data/literature_stats.csv will have, plus the two the review workflow adds.
LIT_COLUMNS = ["lit_id", "entry_id", "scope_type", "scope_id", "effect", "readout",
               "n_papers_total", "n_papers_using", "pmcids", "evidence", "confidence",
               "source", "retrieved_on", "method_version", "notes"]
REVIEW_COLUMNS = ["decision", "your_note"]

# A readout claim resting on one paper is an anecdote. Two is the least that should reach the site,
# and matches readout_gaps.py's own default.
MIN_PAPERS = 2
# A share needs a denominator. With one paper extracted, "1 of 1 papers" is a literally true 100%
# that would feed the weighted prevalence score at weight 0.3 and read as certainty on the site.
# Five is the least denominator worth publishing a percentage from. This is why the structured
# extraction path needs about ten papers per protein before it contributes anything, while the
# readout-gap path works off the whole fetched set from the first run.
MIN_PAPERS_TOTAL = 5


def log(msg):
    print(msg)
    sys.stdout.flush()


def read_csv_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))


def read_json(path, default=None):
    try:
        with open(path, encoding="utf-8") as h:
            return json.load(h)
    except (OSError, ValueError):
        return default


def split(text):
    return [t.strip() for t in re.split(r"[;|]", text or "") if t.strip()]


# ---------------------------------------------------------------- matching papers to playbook entries

def load_entries():
    """Playbook entries with the function labels they answer to and whether they are cell-based.
    Same matching rule build_database.using_for applies, so a row generated here lands on the row
    the build would generate."""
    out = []
    for a in read_csv_rows(os.path.join(DATA, "assays.csv")):
        labels = set(split(a.get("evidence_labels") or a.get("evidence_function")))
        if not labels:
            continue
        scope = (a.get("scope") or "").lower()
        out.append({"entry_id": a["entry_id"], "labels": labels,
                    "classes": set(split(a.get("classes"))),
                    "cell": "cell-based" in scope and "free" not in scope,
                    "readouts": {r.lower() for r in split(a.get("readouts"))}})
    return out


def entries_for(entries, classes, function_labels, cell):
    """Entries whose function labels overlap what the paper measured, in the right context."""
    want = set(function_labels)
    for e in entries:
        if not (e["classes"] & (set(classes) | {"generic"})):
            continue
        if e["cell"] != bool(cell):
            continue
        if e["labels"] & want:
            yield e


def classes_of(gene, survey):
    row = survey.get(gene) or {}
    return [c for c in split(row.get("classes")) if c]


# usage_stats.csv has no "any" effect: build_database writes one row per effect VIEW (inhibition,
# activation, inhibition_labeled, activation_with_unspecified), and a literature row keyed on "any"
# would match nothing and do nothing, silently.
#
# Most paper findings carry no direction. An SPR measurement or a thermal shift is the same
# experiment whichever way the compound acts, so such a finding belongs under every view. Writing
# it out once per view made the review file four times longer for no decision gained, so it is
# proposed once with effect "*" and build_database expands it.
ALL_EFFECTS = "*"


# ---------------------------------------------------------------- reading the extractions

def best_quote(assay, field="technology_labels"):
    """The verbatim sentence behind a claim. `evidence` is a list of {field, quote, section}; the
    one naming the field we are acting on is the one a reviewer needs, so prefer it."""
    ev = assay.get("evidence")
    if isinstance(ev, dict):
        ev = [ev]
    best = ""
    for item in ev or []:
        if not isinstance(item, dict):
            continue
        q = " ".join((item.get("quote") or "").split())
        if not q:
            continue
        if (item.get("field") or "") == field:
            return q
        best = best or q
    return best


def extraction_rows(paths, entries, survey):
    """One record per (gene, entry, effect, readout), counting the papers that support it.

    An assay in a paper contributes its technology_labels to the entries whose function it matches,
    in its own context. This is the same pairing rule the campaign side uses: a binding panel in the
    same paper must not credit its readout to the enzyme-activity entry.

    The gene is taken per assay, from assay.target.gene_symbol, not per paper: a paper routinely
    assays several proteins, and the whole point of the pairing rule is that each readout is
    credited to the protein it measured."""
    seen = collections.defaultdict(lambda: {"papers": set(), "quotes": [], "conf": []})
    totals = collections.defaultdict(set)   # gene -> papers that yielded any usable assay
    for path in paths:
        doc = read_json(path)
        if not doc:
            continue
        paper = doc.get("paper") or {}
        pmcid = paper.get("pmcid") or os.path.splitext(os.path.basename(path))[0]
        fallback = (paper.get("target_genes") or [None])[0]
        for assay in doc.get("assays") or []:
            gene = ((assay.get("target") or {}).get("gene_symbol")) or fallback
            fns = assay.get("function_labels") or []
            techs = assay.get("technology_labels") or []
            if not gene or not fns or not techs:
                continue
            classes = classes_of(gene, survey)
            scope = (assay.get("scope") or "").lower()
            cell = "cell" in scope and "free" not in scope
            effect = (assay.get("effect_direction") or "").strip().lower()
            effect = effect if effect in ("inhibition", "activation") else "any"
            quote = best_quote(assay)
            for e in entries_for(entries, classes, fns, cell):
                totals[gene].add(pmcid)
                for t in techs:
                    rec = seen[(gene, e["entry_id"], effect, t)]
                    rec["papers"].add(pmcid)
                    if quote and len(rec["quotes"]) < 3:
                        rec["quotes"].append("%s: %s" % (pmcid, quote[:160]))
                    rec["conf"].append(assay.get("confidence") or "")
                # the entry-level row: this paper measured this entry's function at all
                seen[(gene, e["entry_id"], effect, "")]["papers"].add(pmcid)
    return seen, {g: len(v) for g, v in totals.items()}


def gap_rows(entries, survey):
    """Readouts readout_gaps.py found in papers that no campaign for the protein states. These are
    the rows the exercise exists for, so they are proposed even when no structured extraction has
    been filled for the paper yet."""
    out = {}
    for r in read_csv_rows(os.path.join(PAPERS, "readout_gaps.csv")):
        gene, readout = r.get("gene"), r.get("readout_new_in_papers")
        if not gene or not readout:
            continue
        try:
            n_using = int(r.get("papers_naming_it") or 0)
            n_total = int(r.get("papers_searched") or 0)
        except ValueError:
            continue
        for eid in split(r.get("entries_that_would_gain")):
            # "" is not an effect view; build_rows expands it to the views the protein has.
            out[(gene, eid, "", readout)] = {"n_total": n_total, "n_using": n_using,
                                             "source": "readout_gaps.py"}
    return out


# ---------------------------------------------------------------- writing the chunk

def build_rows(args):
    survey = {r["gene"]: r for r in read_csv_rows(os.path.join(PAPERS, "survey.csv"))}
    entries = load_entries()
    by_id = {e["entry_id"]: e for e in entries}
    want = {g.strip().upper() for g in (args.genes or "").split(",") if g.strip()}

    paths = sorted(p for p in glob.glob(os.path.join(PAPERS, "PMC*.json")))
    seen, totals = extraction_rows(paths, entries, survey)
    gaps = gap_rows(entries, survey)
    log("  %d structured extraction(s), %d gap row(s)" % (len(paths), len(gaps)))

    rows, skipped = [], collections.Counter()
    keys = sorted(set(seen) | set(gaps))
    for gene, eid, effect, readout in keys:
        if want and gene.upper() not in want:
            continue
        rec = seen.get((gene, eid, effect, readout))
        gap = gaps.get((gene, eid, effect, readout))
        if rec:
            n_using, n_total = len(rec["papers"]), max(totals.get(gene, 0), len(rec["papers"]))
            pmcids, quotes = sorted(rec["papers"]), rec["quotes"]
            src = "paper extraction (schema/paper_extraction.schema.json)"
            conf = min(rec["conf"] or ["medium"], key=lambda c: {"high": 2, "medium": 1}.get(c, 0))
        else:
            n_using, n_total = gap["n_using"], gap["n_total"]
            pmcids, quotes = [], []
            src = "readout_gaps.py (methods text near a subject mention of the protein)"
            conf = "medium"
        if readout and n_using < MIN_PAPERS:
            skipped["readout seen in fewer than %d papers" % MIN_PAPERS] += 1
            continue
        if n_total < MIN_PAPERS_TOTAL:
            skipped["denominator below %d papers" % MIN_PAPERS_TOTAL] += 1
            continue
        if eid not in by_id:
            skipped["unknown_entry"] += 1
            continue
        # A readout the entry already lists is confirmation, not news. Still proposed, because the
        # per-protein share is what the site shows, but flagged so the reviewer can tell them apart.
        already = readout and readout.lower() in by_id[eid]["readouts"]
        # A finding that named a direction keeps it; one that did not applies to every view.
        eff = effect if effect in ("inhibition", "activation") else ALL_EFFECTS
        rows.append({
            "lit_id": "%s-%s-%s%s" % (gene, eid, "all" if eff == ALL_EFFECTS else eff[:3],
                                      "-" + re.sub(r"\W+", "", readout)[:12] if readout else ""),
            "entry_id": eid, "scope_type": "gene", "scope_id": gene, "effect": eff,
            "readout": readout, "n_papers_total": n_total, "n_papers_using": n_using,
            "pmcids": "; ".join(pmcids), "evidence": " | ".join(quotes),
            "confidence": conf, "source": src, "retrieved_on": args.today,
            "method_version": "lit-" + args.today.replace("-", ""),
            "notes": ("readout already listed on this entry; adds the protein-specific share"
                      if already else
                      "readout NOT listed on this entry for this protein" if readout else
                      "entry-level: papers measuring this entry's function")
                     + ("; no direction in the papers, so it counts under every effect view"
                        if eff == ALL_EFFECTS else ""),
            "decision": "", "your_note": ""})
    return rows, skipped


def write_chunk(rows, args):
    os.makedirs(args.chunk, exist_ok=True)
    path = os.path.join(args.chunk, "literature_proposals.csv")
    with open(path, "w", newline="", encoding="utf-8-sig") as h:
        w = csv.DictWriter(h, fieldnames=LIT_COLUMNS + REVIEW_COLUMNS, lineterminator="\r\n")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["scope_id"], r["entry_id"], r["readout"])):
            w.writerow(r)
    by_gene = collections.Counter(r["scope_id"] for r in rows)
    new = [r for r in rows if r["readout"] and "NOT listed" in r["notes"]]
    md = os.path.join(args.chunk, "REVIEW.md")
    with open(md, "w", encoding="utf-8") as h:
        h.write(CHUNK_README.replace("__N__", str(len(rows)))
                .replace("__NGENES__", str(len(by_gene)))
                .replace("__NNEW__", str(len(new)))
                .replace("__MINP__", str(MIN_PAPERS))
                .replace("__TABLE__", summary_table(rows)))
    return path, md


def summary_table(rows):
    out = ["| protein | entry | effect | readout | papers | of | what it adds |",
           "|---|---|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda r: (r["scope_id"], r["entry_id"], r["effect"], r["readout"]))[:60]:
        adds = ("NEW readout for this protein" if r["readout"] and "NOT listed" in r["notes"]
                else "protein-specific share" if r["readout"] else "function occurrence")
        out.append("| %s | %s | %s | %s | %s | %s | %s |"
                   % (r["scope_id"], r["entry_id"], r["effect"], r["readout"] or "(whole entry)",
                      r["n_papers_using"], r["n_papers_total"], adds))
    return "\n".join(out)


CHUNK_README = """# Literature evidence from mined papers

__N__ proposed row(s) across __NGENES__ protein(s). __NNEW__ of them are readouts that papers show
and **no campaign for that protein states** - the case this exercise exists for.

## What accepting these does

They become `data/literature_stats.csv`, which `scripts/build_database.py` reads while generating
`data/usage_stats.csv`. It fills the two columns that have always been there and always been empty,
`n_literature_total` and `n_literature_using`. Then:

- `tools/build_db.py` already scores them: the prevalence score is a weighted mean of the campaign
  share and the literature share, weights in `data/config.json` (campaigns 0.7, literature 0.3).
- `docs/index.html` already renders them: the Function occurrence line gains "; 60% of 25 papers".
- A readout with no campaign behind it appears with a campaign count of zero and a paper count that
  is not zero, which is an honest statement of what is known.

Nothing in the database schema or the site had to change to accept this. The columns, the scorer and
the display were all built for literature evidence and have been waiting for data.

## What it does NOT do

It does not touch campaign counts. `n_campaigns_using` keeps coming from ChEMBL, PubChem and
BindingDB on every rebuild. Paper evidence is a second, separately weighted channel, never a
correction to the first. That is why these rows live in their own file and are folded in at
generation time, instead of being written into `usage_stats.csv` by hand: a hand-written row there
is preserved by the merge and would freeze that protein's campaign counts permanently.

## How to review

Set `decision` on each row to `accept`, `reject`, or `edit` (with the corrected value in the row and
a note). Every row carries `pmcids` and a verbatim `evidence` quote, so a claim can be checked
against the paper without running anything.

Rows resting on fewer than __MINP__ papers were dropped before this file was written.

Then:

    py tools\\apply_review.py __CHUNK__ --apply
    py scripts\\build_database.py --aggregate-only --merge
    py tools\\build_db.py

## Proposed rows

__TABLE__
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chunk", default=os.path.join("review", "08_literature"),
                    help="review chunk folder to write (default review\\08_literature)")
    ap.add_argument("--genes", default="", help="comma-separated subset")
    ap.add_argument("--today", default=None, help="date stamp (default: today)")
    ap.add_argument("--dry-run", action="store_true", help="print, write nothing")
    args = ap.parse_args(argv)
    if not args.today:
        from datetime import date
        args.today = date.today().isoformat()

    log("Building literature proposals...")
    rows, skipped = build_rows(args)
    log("  %d proposed row(s)" % len(rows))
    for k, n in skipped.most_common():
        log("  %d skipped (%s)" % (n, k))
    if not rows:
        log("  nothing to propose")
        return 0
    if args.dry_run:
        log("\n" + summary_table(rows))
        return 0
    path, md = write_chunk(rows, args)
    with open(md, encoding="utf-8") as h:
        text = h.read().replace("__CHUNK__", args.chunk)
    with open(md, "w", encoding="utf-8") as h:
        h.write(text)
    log("  %s" % os.path.relpath(path, ROOT))
    log("  %s" % os.path.relpath(md, ROOT))
    log("\nReview the decision column, then:")
    log("  py tools\\apply_review.py %s --apply" % args.chunk)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
