#!/usr/bin/env python
"""Aggregate the crawl's paper records into candidate changes to the curated playbook.

    py tools\\papers_aggregate.py                  (all batches)
    py tools\\papers_aggregate.py --batches 2,3    (only these)
    py tools\\papers_aggregate.py --min-papers 3   (stricter threshold)

Why this exists
---------------
scripts/literature_stats.py already turns paper records into data/literature_stats.csv: the
*counting* tier, "60% of 25 papers", per protein. It cannot touch data/assays.csv, so none of the
findings that are about the playbook itself -- a readout no entry lists, a function no entry
covers -- can reach the app through it. This file is that other half.

The governing rule, from schema/PAPER_INTEGRATION.md: entry-level columns are per FUNCTION CLASS,
not per protein. data/assays.csv has ~110 entries keyed by 37 classes, so a KDR paper does not
improve KDR, it improves `protein_kinase` and every kinase inherits it. Aggregating per protein
would propose the same change 18 times and make a kinase-heavy batch look like 18 independent
observations. So every candidate here is keyed by (class, ...) and carries how many distinct
proteins and papers stand behind it.

What it reads
-------------
crawl/evidence.csv      one row per accepted assay record, all batches, with protein_classes,
                        function_labels, technology_labels, scope and the matched entry ids
                        already resolved by scripts/paper_crawl.py log.
crawl/term_catalog.csv  the vocabulary candidates, with status candidate/supported.
data/assays.csv         what the playbook holds today, to subtract.

What it writes, into build_output/papers/
-----------------------------------------
readout_candidates.csv     (class, entry_id, readout) a technology no entry of that class lists
new_entry_candidates.csv   (class, function, scope) records that matched NO entry at all -- the
                           playbook has no row for this, which is the expensive kind of gap
term_candidates.csv        vocabulary terms at or above the threshold, from the term catalog
link_candidates.csv        (gene, entry_id) informational: entries a protein is not linked to
                           today. Per protein by nature, and literature_stats.py already carries
                           this to the site, so it is a report here, not a proposal.

Nothing here writes to data/ or to review/. tools/papers_to_review.py renders a chunk from these.

Thresholds
----------
A candidate needs --min-papers distinct papers (default 2: one paper is an anecdote). Records with
confidence `low` are excluded from counting and reported separately, per PAPER_INTEGRATION.md.
"""
import argparse
import collections
import csv
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
CRAWL = os.path.join(ROOT, "crawl")
OUT = os.path.join(ROOT, "build_output", "papers")

sys.path.insert(0, os.path.join(ROOT, "scripts"))
try:
    import readout_gaps as rg
    GENERIC = {g.lower() for g in rg.GENERIC_TECHNIQUES}
except Exception:                                    # keep working if that module moves
    GENERIC = {"western blot", "mass spectrometry", "qpcr", "elisa / immunosorbent",
               "flow cytometry", "microscopy / imaging"}

MIN_PAPERS = 2

# ---------------------------------------------------------------- readout identity
#
# A readouts cell is curated free text ("HTRF / TR-FRET cAMP"), while a proposal is a controlled
# label ("TR-FRET / HTRF"). Comparing them literally proposed 7 readouts in 27 that the entry
# already listed, measured on review/11_papers_batches_0_4. Word order, punctuation, a trailing
# example and a house abbreviation all defeated the match. So both sides are reduced to a set of
# meaningful tokens and compared by subset, in either direction: a proposal that says no more than
# the cell already says is not a gap.

# Expanded before tokenising, longest first, on word boundaries. Only abbreviations that actually
# collided; this is a lookup table, not a theory of assay naming.
READOUT_ALIASES = [
    ("kinase-glo", "luminescent atp depletion"),     # Kinase-Glo reads the ATP left unconsumed
    ("nanobret", "bret"),                            # a BRET implementation, not a separate readout
    ("htrf", "tr fret"),                             # homogeneous time-resolved fluorescence IS TR-FRET
    ("fp", "fluorescence polarization"),
    ("dsf", "thermal shift"),
    ("bli", "biolayer interferometry"),
]

# Dropped as carrying no discriminating meaning. "binding" and "detection" matter here: they are
# what separates "Radioligand binding" from "Radioligand filtration or scintillation proximity".
# `fluorescence` and `luminescence` are deliberately NOT dropped: without them "Fluorescence
# polarization" collapses to one token, and one-token keys match far too much.
READOUT_STOPWORDS = {
    "a", "an", "and", "or", "the", "of", "in", "for", "with", "to", "by", "on", "e", "g", "i",
    "assay", "assays", "based", "format", "readout", "readouts", "detection", "detect", "binding",
    "label", "labeled", "labelled", "free", "competition", "competitive", "cell", "cells",
}

# A one-token key is usually too weak to prove identity: `luciferase` alone would otherwise match
# `split-luciferase complementation`, which is a different assay. A single token may stand for a
# whole readout only if it is specific enough on its own.
SPECIFIC_TOKENS = {
    "radioligand", "transcreener", "lanthascreen", "kinomescan", "itc", "spr", "cetsa",
    "alphascreen", "alphalisa", "crystallography", "electrophysiology", "thermophoresis",
}


def readout_key(text):
    """A readout reduced to the tokens that identify it. Returns a frozenset."""
    s = (text or "").lower()
    s = re.sub(r"\(.*?\)", " ", s)                   # trailing examples are not part of the identity
    for short, long in READOUT_ALIASES:
        s = re.sub(r"\b%s\b" % re.escape(short), long, s)
    toks = [t for t in re.split(r"[^a-z0-9]+", s) if t]
    return frozenset(t for t in toks if t not in READOUT_STOPWORDS)


def _same_readout(a, b):
    """Identical, or one says no more than the other and is specific enough to be trusted."""
    if not a or not b:
        return False
    if a == b:
        return True
    small, large = (a, b) if len(a) < len(b) else (b, a)
    if not small <= large:
        return False
    return len(small) >= 2 or next(iter(small)) in SPECIFIC_TOKENS


def already_listed(tech, listed_keys):
    """-> the readout the cell already lists that says the same thing, or None."""
    k = readout_key(tech)
    if not k:
        return None
    for original, lk in listed_keys:
        if _same_readout(k, lk):
            return original
    return None


def log(msg):
    print(msg)
    sys.stdout.flush()


def read_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))


def split(text):
    return [t.strip() for t in re.split(r"[;|]", text or "") if t.strip()]


def write_rows(path, cols, rows):
    with open(path, "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


# ---------------------------------------------------------------- the playbook as it stands

def load_playbook():
    """-> (entries by id, readouts already listed per class, function labels covered per class).

    `covered` answers "does any entry of this class, in this context, already measure this
    function" -- the test that decides whether a record is a readout gap (entry exists, readout
    missing) or an entry gap (no row at all). Context matters: a cell-based record must not be
    answered by a cell-free entry, which is the same pairing rule literature_stats.entries_for
    and build_database.using_for apply.
    """
    entries, by_class_readouts, covered = {}, collections.defaultdict(dict), collections.defaultdict(set)
    controlled = set()
    for a in read_rows(os.path.join(DATA, "assays.csv")):
        scope = (a.get("scope") or "").lower()
        cell = "cell-based" in scope and "free" not in scope
        e = {"entry_id": a["entry_id"], "name": a.get("assay_name") or "",
             "classes": split(a.get("classes")), "cell": cell,
             "category": a.get("assay_category") or "",
             "readouts": split(a.get("readouts")),
             "labels": set(split(a.get("evidence_labels") or a.get("evidence_function")))}
        entries[e["entry_id"]] = e
        controlled.update(r.lower() for r in e["readouts"])
        controlled.update(lab.lower() for lab in e["labels"])
        for c in e["classes"]:
            for r in e["readouts"]:
                # original text -> its identity key, so a suppression can name what it matched
                by_class_readouts[(c, cell)][r] = readout_key(r)
            for lab in e["labels"]:
                covered[(c, cell)].add(lab.lower())
    return entries, by_class_readouts, covered, controlled


# ---------------------------------------------------------------- aggregation

def aggregate(ev, entries, by_class_readouts, covered, min_papers):
    """Walk the records once, bucketing each into the gap it represents.

    One record can be several things at once: it can name a readout the class does not list AND a
    function no entry covers. Both are recorded, because they are different proposals with
    different costs -- a readout is one cell in an existing row, a new entry is a new row the owner
    has to rate.
    """
    readout = collections.defaultdict(lambda: {"papers": set(), "genes": set(), "quote": "", "entries": set()})
    newentry = collections.defaultdict(lambda: {"papers": set(), "genes": set(), "quote": "",
                                                "techs": collections.Counter()})
    link = collections.defaultdict(lambda: {"papers": set(), "quote": ""})
    incoherent = {}
    suppressed = collections.Counter()      # (class, scope, proposed, matched) -> n records
    dropped_low, no_class = 0, set()

    for r in ev:
        if (r.get("confidence") or "").strip().lower() == "low":
            dropped_low += 1
            continue
        gene, pm = r.get("gene") or "", r.get("pmcid") or ""
        classes = [c for c in split(r.get("protein_classes")) if c and c != "none"]
        if not classes:
            no_class.add(gene)
            continue
        cell = (r.get("scope") or "").strip() == "cell-based"
        quote = (r.get("quote") or "")[:300]
        matched = split(r.get("entries"))

        # 1. readouts: a technology this class does not list anywhere, in this context
        for tech in split(r.get("technology_labels")):
            if tech.lower() in GENERIC:
                continue
            for c in classes:
                listed = by_class_readouts.get((c, cell), {})
                hit = already_listed(tech, listed.items())
                if hit is not None:
                    # not a gap: the cell says this already, in different words
                    suppressed[(c, "cell-based" if cell else "cell-free", tech, hit)] += 1
                    continue
                k = (c, "cell-based" if cell else "cell-free", tech)
                d = readout[k]
                d["papers"].add(pm)
                d["genes"].add(gene)
                d["entries"].update(matched)
                d["quote"] = d["quote"] or quote

        # 2. no entry matched at all: the playbook has no row for this function in this context
        if not matched:
            for c in classes:
                for lab in split(r.get("function_labels")):
                    if lab.lower() in covered.get((c, cell), set()):
                        continue          # a row exists; the record failed to match for another reason
                    # A function that says "cellular" cannot be a cell-free gap: the record
                    # contradicts itself, so it is a flaw in the extraction, not a hole in the
                    # playbook. Proposing an entry for it would put that contradiction in the
                    # database. These come from lysate methods (DARTS, CETSA on lysate) recorded as
                    # cell-free while keeping the cellular-engagement function label.
                    if not cell and "cellular" in lab.lower():
                        incoherent[(lab, gene, pm)] = (r.get("technology_labels") or "")
                        continue
                    k = (c, "cell-based" if cell else "cell-free", lab)
                    d = newentry[k]
                    d["papers"].add(pm)
                    d["genes"].add(gene)
                    d["quote"] = d["quote"] or quote
                    for tech in split(r.get("technology_labels")):
                        d["techs"][tech] += 1

        # 3. informational: entries this protein is not shown for today
        for e in split(r.get("new_entry_links")):
            d = link[(gene, e)]
            d["papers"].add(pm)
            d["quote"] = d["quote"] or quote

    rows_readout = [{"class_id": c, "scope": sc, "readout": t,
                     "n_papers": len(d["papers"]), "n_proteins": len(d["genes"]),
                     "entries_seen": "; ".join(sorted(d["entries"])),
                     "proteins": "; ".join(sorted(d["genes"])[:8]),
                     "pmcids": "; ".join(sorted(d["papers"])[:10]), "quote": d["quote"]}
                    for (c, sc, t), d in readout.items() if len(d["papers"]) >= min_papers]
    rows_entry = [{"class_id": c, "scope": sc, "function_label": lab,
                   "n_papers": len(d["papers"]), "n_proteins": len(d["genes"]),
                   "readouts_seen": "; ".join(t for t, _ in d["techs"].most_common(6)),
                   "proteins": "; ".join(sorted(d["genes"])[:8]),
                   "pmcids": "; ".join(sorted(d["papers"])[:10]), "quote": d["quote"]}
                  for (c, sc, lab), d in newentry.items() if len(d["papers"]) >= min_papers]
    rows_link = [{"gene": g, "entry_id": e, "n_papers": len(d["papers"]),
                  "pmcids": "; ".join(sorted(d["papers"])[:10]), "quote": d["quote"]}
                 for (g, e), d in link.items() if len(d["papers"]) >= min_papers]
    rows_bad = [{"function_label": lab, "gene": g, "pmcid": pm, "scope": "cell-free",
                 "technology_labels": tech,
                 "problem": "function says cellular but the record is cell-free"}
                for (lab, g, pm), tech in sorted(incoherent.items())]
    rows_sup = [{"class_id": c, "scope": sc, "proposed_readout": t, "already_listed_as": hit,
                 "n_records": n}
                for (c, sc, t, hit), n in sorted(suppressed.items(), key=lambda kv: -kv[1])]
    return rows_readout, rows_entry, rows_link, rows_bad, rows_sup, dropped_low, no_class


def controlled_vocabulary():
    """Everything the database already knows a term by, from both places that define it.

    scripts/term_table.csv is the authority: it is what scripts/m2_label.py matches campaign text
    against, so it is the table a newly accepted term would join, and its `synonyms` column means a
    term can already be covered under another name ("thallium flux" is a synonym of the function
    "Ion channel activity"). data/assays.csv is checked as well because a value can be in use as a
    readout on an entry without being a labeller term.
    """
    out = set()
    for t in read_rows(os.path.join(ROOT, "scripts", "term_table.csv")):
        out.add((t.get("term") or "").strip().lower())
        out.update(s.strip().lower() for s in split(t.get("synonyms")))
    for a in read_rows(os.path.join(DATA, "assays.csv")):
        out.update(r.lower() for r in split(a.get("readouts")))
        out.update(l.lower() for l in split(a.get("evidence_labels") or a.get("evidence_function")))
    return {x for x in out if x}


def aggregate_systems(ev, min_papers):
    """Candidate values for the `system` column: the biological preparation an assay ran in.

    These do not come from a record field -- the crawl's extraction schema has no slot for the
    preparation, which is exactly the gap. They come from the term catalog, where the crawl filed
    them under the kinds `system` and `counter_screen`, and are attributed back to a class and an
    entry by joining each occurrence to its record on assay_id.

    Expect these to sit BELOW the threshold for a long time. A preparation is named once per paper
    and papers rarely repeat each other's, so corroboration needs several papers in the same class
    rather than more papers overall -- an argument for choosing batches by empty class.
    """
    kinds = {}
    for t in read_rows(os.path.join(CRAWL, "term_catalog.csv")):
        if (t.get("kind") or "") in ("system", "counter_screen"):
            kinds[(t.get("proposed_term") or "").strip()] = t.get("kind")
    if not kinds:
        return [], []
    by_record = {r.get("assay_id"): r for r in ev}
    seen = collections.defaultdict(lambda: {"papers": set(), "genes": set(), "quote": "",
                                            "entries": set()})
    for o in read_rows(os.path.join(CRAWL, "term_occurrences.csv")):
        term = (o.get("catalog_term") or "").strip()
        if term not in kinds:
            continue
        rec = by_record.get(o.get("assay_id"))
        if not rec:
            continue            # the record was skipped or is below confidence; nothing to attach to
        classes = [c for c in split(rec.get("protein_classes")) if c and c != "none"]
        scope = "cell-based" if (rec.get("scope") or "").strip() == "cell-based" else "cell-free"
        for c in classes:
            d = seen[(c, scope, term, kinds[term])]
            d["papers"].add(o.get("pmcid") or "")
            d["genes"].add(o.get("gene") or "")
            d["entries"].update(split(rec.get("entries")))
            d["quote"] = d["quote"] or (o.get("quote") or "")[:300]
    rows = [{"class_id": c, "scope": sc, "system": term, "kind": kind,
             "n_papers": len(d["papers"]), "n_proteins": len(d["genes"]),
             "entries_seen": "; ".join(sorted(d["entries"])),
             "proteins": "; ".join(sorted(x for x in d["genes"] if x)),
             "pmcids": "; ".join(sorted(x for x in d["papers"] if x)), "quote": d["quote"]}
            for (c, sc, term, kind), d in seen.items()]
    ready = [r for r in rows if r["n_papers"] >= min_papers]
    waiting = [r for r in rows if r["n_papers"] < min_papers]
    return ready, waiting


def aggregate_terms(min_papers, controlled):
    """Vocabulary terms worth proposing, with how many papers actually used each.

    The catalog's own `status` column says candidate (1 paper) or supported (2+), but it was set by
    hand as the crawl went, so the occurrence table is the authority on the count.

    A term that is ALREADY a controlled value is dropped. The catalog accumulated a few of these by
    being edited as the crawl ran: T005 began as "imaging flow cytometry", was folded into the
    existing "Flow cytometry" value, and kept the catalog row -- so it would be proposed as new
    vocabulary when the database has held it all along. Proposing an existing term is worse than
    noise: accepting it implies an m2 relabel that would change nothing.
    """
    occ = read_rows(os.path.join(CRAWL, "term_occurrences.csv"))
    seen = collections.defaultdict(lambda: {"papers": set(), "genes": set(), "quote": ""})
    for o in occ:
        d = seen[o.get("catalog_term") or ""]
        d["papers"].add(o.get("pmcid") or "")
        d["genes"].add(o.get("gene") or "")
        d["quote"] = d["quote"] or (o.get("quote") or "")[:300]
    out, already = [], []
    for t in read_rows(os.path.join(CRAWL, "term_catalog.csv")):
        name = (t.get("proposed_term") or "").strip()
        d = seen.get(name, {"papers": set(), "genes": set(), "quote": ""})
        if len(d["papers"]) < min_papers:
            continue
        if name.lower() in controlled:
            already.append(name)
            continue
        out.append({"term_id": t.get("term_id"), "proposed_term": t.get("proposed_term"),
                    "kind": t.get("kind"), "nearest_controlled_value": t.get("nearest_controlled_value"),
                    "definition": t.get("definition"), "variants": t.get("variants"),
                    # carried straight through to scripts/term_table.csv by papers_to_review.py
                    "label_synonyms": t.get("label_synonyms") or "",
                    "label_suppresses": t.get("label_suppresses") or "",
                    "n_papers": len(d["papers"]), "n_proteins": len(d["genes"]),
                    "proteins": "; ".join(sorted(x for x in d["genes"] if x)[:8]),
                    "pmcids": "; ".join(sorted(x for x in d["papers"] if x)[:10]),
                    "quote": d["quote"], "catalog_status": t.get("status")})
    return out, already


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--batches", help="comma-separated batch numbers (default: all)")
    ap.add_argument("--min-papers", type=int, default=MIN_PAPERS,
                    help="distinct papers a candidate needs (default %d)" % MIN_PAPERS)
    args = ap.parse_args(argv)

    ev = read_rows(os.path.join(CRAWL, "evidence.csv"))
    if not ev:
        raise SystemExit("crawl/evidence.csv is empty or missing: run paper_crawl.py log first")
    if args.batches:
        want = {b.strip() for b in args.batches.split(",") if b.strip()}
        ev = [r for r in ev if (r.get("batch") or "").strip() in want]
        if not ev:
            raise SystemExit("no records for batch(es) %s" % args.batches)

    entries, by_class_readouts, covered, controlled = load_playbook()
    log("Aggregating %d record(s) against %d playbook entries, threshold %d paper(s)."
        % (len(ev), len(entries), args.min_papers))

    ro, ne, lk, bad, sup, dropped_low, no_class = aggregate(ev, entries, by_class_readouts, covered,
                                                            args.min_papers)
    terms, already = aggregate_terms(args.min_papers, controlled | controlled_vocabulary())
    systems, systems_waiting = aggregate_systems(ev, args.min_papers)

    os.makedirs(OUT, exist_ok=True)
    ro.sort(key=lambda r: (-r["n_papers"], r["class_id"], r["readout"]))
    ne.sort(key=lambda r: (-r["n_papers"], r["class_id"], r["function_label"]))
    lk.sort(key=lambda r: (-r["n_papers"], r["gene"]))
    terms.sort(key=lambda r: (-r["n_papers"], r["proposed_term"]))

    write_rows(os.path.join(OUT, "readout_candidates.csv"),
               ["class_id", "scope", "readout", "n_papers", "n_proteins", "entries_seen",
                "proteins", "pmcids", "quote"], ro)
    write_rows(os.path.join(OUT, "new_entry_candidates.csv"),
               ["class_id", "scope", "function_label", "n_papers", "n_proteins", "readouts_seen",
                "proteins", "pmcids", "quote"], ne)
    write_rows(os.path.join(OUT, "term_candidates.csv"),
               ["term_id", "proposed_term", "kind", "nearest_controlled_value", "definition",
                "variants", "label_synonyms", "label_suppresses", "n_papers", "n_proteins",
                "proteins", "pmcids", "quote", "catalog_status"], terms)
    write_rows(os.path.join(OUT, "link_candidates.csv"),
               ["gene", "entry_id", "n_papers", "pmcids", "quote"], lk)
    write_rows(os.path.join(OUT, "record_problems.csv"),
               ["function_label", "gene", "pmcid", "scope", "technology_labels", "problem"], bad)
    write_rows(os.path.join(OUT, "readout_suppressed.csv"),
               ["class_id", "scope", "proposed_readout", "already_listed_as", "n_records"], sup)
    syscols = ["class_id", "scope", "system", "kind", "n_papers", "n_proteins", "entries_seen",
               "proteins", "pmcids", "quote"]
    write_rows(os.path.join(OUT, "system_candidates.csv"), syscols,
               sorted(systems, key=lambda r: (-r["n_papers"], r["class_id"])))
    write_rows(os.path.join(OUT, "system_candidates_below_threshold.csv"), syscols,
               sorted(systems_waiting, key=lambda r: (r["class_id"], r["system"])))

    log("")
    log("  readout_candidates.csv     %4d  (class lists no such readout)" % len(ro))
    log("  new_entry_candidates.csv   %4d  (no playbook row covers this function)" % len(ne))
    log("  term_candidates.csv        %4d  (vocabulary at or above threshold)" % len(terms))
    log("  link_candidates.csv        %4d  (informational, per protein)" % len(lk))
    log("  record_problems.csv        %4d  (self-contradictory records, not proposed)" % len(bad))
    if sup:
        log("  readout_suppressed.csv     %4d  (the cell already says this in different words)" % len(sup))
    log("  system_candidates.csv     %4d  (biological preparation, at threshold)" % len(systems))
    if systems_waiting:
        log("  %d system/counter-screen value(s) below the threshold, each named by only one paper;"
            % len(systems_waiting))
        log("  listed in system_candidates_below_threshold.csv. A preparation is stated once per")
        log("  paper, so these corroborate only as more papers of the SAME CLASS arrive.")
    if already:
        log("  %d term(s) dropped: already a controlled value (%s)."
            % (len(already), ", ".join(sorted(set(already))[:5])))
    if dropped_low:
        log("  %d record(s) excluded: confidence low." % dropped_low)
    if no_class:
        log("  %d record(s) carry no protein class, so no entry can receive them: %s"
            % (len(no_class), ", ".join(sorted(no_class))))
    log("")
    log("Next: py tools\\papers_to_review.py --chunk review\\NN_papers_<topic>")
    return 0


if __name__ == "__main__":
    sys.exit(main())
