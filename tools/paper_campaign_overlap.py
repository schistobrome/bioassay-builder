#!/usr/bin/env python
"""Does a mined paper already count as one of the campaigns?

    py scripts\\m2b_fetch_abstracts.py --crawled      build the map first (cached, run once)
    py tools\\paper_campaign_overlap.py

Why this matters
----------------
schema/PAPER_INTEGRATION.md calls this the first thing that will go wrong if ignored. A campaign is
keyed `doc:CHEMBL1155814`; a mined paper carries a PMID and a DOI. If a paper we mined IS one of the
ChEMBL documents behind a campaign, then counting it as literature evidence counts the same
experiment twice -- and until this ran there was no map on disk to tell.

The separate-tier design keeps the current numbers honest either way: papers never enter the
campaign denominator, so no double counting can reach the site today. What this answers is whether
the two tiers could ever be combined, and how badly the overlap would bite if they were.

What it reads
-------------
build_cache/m1_data/chembl_docs/*.json   PMID and DOI per ChEMBL document (from m2b)
build_cache/m1_data/chembl_*_assays.csv  which documents belong to which gene's campaigns
build_output/papers/*.json               the mined papers, with their PMID and DOI

What it writes
--------------
build_output/papers/campaign_overlap.csv   one row per mined paper, matched or not

Reading the result: a `same_gene` match is a real double count waiting to happen. An
`other_gene` match is the same paper behind a campaign for a DIFFERENT protein, which is not a
double count for the protein we mined -- papers routinely assay several targets.

Measured 2026-11-18 against 220 mined papers: 1 same_gene match (PMC11215726, JAK2, CHEMBL5500392).
Rare because ChEMBL curates medicinal-chemistry SAR journals (J Med Chem, Bioorg Med Chem Lett) while
the crawl draws on the PMC open-access subset (Nature Communications, PLoS One) -- but the one match
IS a J Med Chem paper that reached PMC via open access, so the corpora are only mostly disjoint.

One caution from finding it: at 9,200 cached documents this printed 0, and the next 2,417 documents
produced the match. A zero from a partial document cache is not an absence of overlap. Check the
cached-document count against the total m2b reports before believing a zero, and re-run after any
batch that changes where papers come from.

Scope: with `--crawled` the cache covers the crawled genes' campaigns, which makes `same_gene`
complete (mined papers are mined for crawled genes) but leaves `other_gene` undercounted -- a paper
behind a campaign for an uncrawled gene shows as no_match. `--all` closes that, and only matters if
the tiers are ever combined.
"""
import csv
import glob
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M1 = os.path.join(ROOT, "build_cache", "m1_data")
DOCS = os.path.join(M1, "chembl_docs")
PAPERS = os.path.join(ROOT, "build_output", "papers")
OUT = os.path.join(PAPERS, "campaign_overlap.csv")


def log(msg):
    print(msg)
    sys.stdout.flush()


def norm_doi(d):
    d = (d or "").strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d


def load_doc_map():
    """-> (pmid -> {doc_id}, doi -> {doc_id}, how many docs carry each identifier)."""
    by_pmid, by_doi, n, with_pmid, with_doi = {}, {}, 0, 0, 0
    for p in glob.glob(os.path.join(DOCS, "*.json")):
        try:
            with open(p, encoding="utf-8") as h:
                d = json.load(h)
        except (OSError, ValueError):
            continue
        n += 1
        did = d.get("document_chembl_id") or os.path.basename(p)[:-5]
        pmid = str(d.get("pubmed_id") or "").strip()
        doi = norm_doi(d.get("doi"))
        if pmid:
            by_pmid.setdefault(pmid, set()).add(did)
            with_pmid += 1
        if doi:
            by_doi.setdefault(doi, set()).add(did)
            with_doi += 1
    return by_pmid, by_doi, n, with_pmid, with_doi


def load_gene_docs():
    """-> doc_id -> {gene}: which gene's campaigns each ChEMBL document sits behind."""
    out = {}
    for p in glob.glob(os.path.join(M1, "chembl_*_assays.csv")):
        gene = os.path.basename(p)[len("chembl_"):-len("_assays.csv")].upper()
        try:
            with open(p, newline="", encoding="utf-8-sig") as h:
                for r in csv.DictReader(h):
                    d = (r.get("document") or "").strip()
                    if d.startswith("CHEMBL"):
                        out.setdefault(d, set()).add(gene)
        except (OSError, ValueError):
            continue
    return out


def load_mined():
    """-> list of mined papers with the gene each was mined for."""
    out = []
    for p in sorted(glob.glob(os.path.join(PAPERS, "*.json"))):
        if os.path.basename(p) == "index.json":
            continue
        try:
            with open(p, encoding="utf-8") as h:
                d = json.load(h)
        except (OSError, ValueError):
            continue
        paper = d.get("paper") or {}
        genes = {str(g).upper() for g in (paper.get("target_genes") or []) if g}
        # The gene actually assayed can differ from the gene the paper was fetched for, which is
        # the point of recording target.gene_symbol per assay; include both.
        for a in d.get("assays") or []:
            g = ((a.get("target") or {}).get("gene_symbol") or "").upper()
            if g:
                genes.add(g)
        out.append({"pmcid": paper.get("pmcid") or os.path.basename(p)[:-5],
                    "pmid": str(paper.get("pmid") or "").strip(),
                    "doi": norm_doi(paper.get("doi")),
                    "year": paper.get("year") or "",
                    "genes": genes})
    return out


def main():
    if not os.path.isdir(DOCS):
        raise SystemExit("no %s\nRun: py scripts\\m2b_fetch_abstracts.py --crawled" % DOCS)
    by_pmid, by_doi, n_docs, with_pmid, with_doi = load_doc_map()
    if not n_docs:
        raise SystemExit("%s is empty: run m2b_fetch_abstracts.py first" % DOCS)
    gene_docs = load_gene_docs()
    mined = load_mined()

    log("ChEMBL documents cached   : %d" % n_docs)
    log("  with a PubMed id        : %d (%.0f%%)" % (with_pmid, 100.0 * with_pmid / n_docs))
    log("  with a DOI              : %d (%.0f%%)" % (with_doi, 100.0 * with_doi / n_docs))
    log("Mined papers              : %d" % len(mined))
    # Without this, a zero result is ambiguous: a join that found nothing looks exactly like a join
    # that had nothing to work with. Print the mined side's identifier coverage too, so "no overlap"
    # can be told apart from "no identifiers".
    m_pmid = sum(1 for m in mined if m["pmid"])
    m_doi = sum(1 for m in mined if m["doi"])
    m_neither = sum(1 for m in mined if not m["pmid"] and not m["doi"])
    log("  with a PubMed id        : %d" % m_pmid)
    log("  with a DOI              : %d" % m_doi)
    if m_neither:
        log("  with NEITHER            : %d  (these cannot match, whatever the truth is)" % m_neither)

    rows, same, other, none = [], 0, 0, 0
    for m in mined:
        hits = set()
        if m["pmid"]:
            hits |= by_pmid.get(m["pmid"], set())
        if m["doi"]:
            hits |= by_doi.get(m["doi"], set())
        hit_genes = set()
        for d in hits:
            hit_genes |= gene_docs.get(d, set())
        overlap = m["genes"] & hit_genes
        verdict = "same_gene" if overlap else ("other_gene" if hits else "no_match")
        same += verdict == "same_gene"
        other += verdict == "other_gene"
        none += verdict == "no_match"
        rows.append({"pmcid": m["pmcid"], "pmid": m["pmid"], "doi": m["doi"], "year": m["year"],
                     "mined_for": "; ".join(sorted(m["genes"])), "verdict": verdict,
                     "chembl_documents": "; ".join(sorted(hits)),
                     "campaign_genes": "; ".join(sorted(hit_genes)),
                     "genes_in_common": "; ".join(sorted(overlap))})

    with open(OUT, "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=["pmcid", "pmid", "doi", "year", "mined_for", "verdict",
                                          "chembl_documents", "campaign_genes", "genes_in_common"])
        w.writeheader()
        w.writerows(rows)

    log("")
    log("  same_gene   %4d  the paper is already a campaign for a protein we mined it for" % same)
    log("  other_gene  %4d  it is a campaign, but for a different protein" % other)
    log("  no_match    %4d  not among the cached ChEMBL documents" % none)
    log("")
    log("Wrote %s" % os.path.relpath(OUT, ROOT))
    if same:
        log("")
        log("Double counting is possible for these; keep the tiers separate:")
        for r in rows:
            if r["verdict"] == "same_gene":
                log("   %-14s %-10s %s" % (r["pmcid"], r["genes_in_common"], r["chembl_documents"]))
    elif m_pmid or m_doi:
        log("No mined paper is a campaign for the same protein, in the documents cached so far.")
        # Learned the hard way: at 9,200 of 12,867 documents this printed 0, and the next 2,417
        # produced a match. Identifiers on both sides prove the join can work; they say nothing
        # about whether the document cache is complete enough for a zero to mean anything.
        log("Identifiers are present on both sides, so the join works -- but a zero only means")
        log("'not among the %d documents cached'. Re-run m2b until it reports no failures." % n_docs)
    else:
        log("INCONCLUSIVE: no mined paper carries a PubMed id or a DOI, so nothing could match.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
