#!/usr/bin/env python3
"""Search Europe PMC (free, no key) for papers about an assay, to find citations and rough counts.

    python scripts/find_literature.py "HDAC1 inhibitor screening fluorogenic" --max 25
    python scripts/find_literature.py "BRD4 AlphaScreen" --count-only

Writes data/literature_candidates.csv (doi, pmid, title, year, journal, cited_by, open_access) for you to read and
pick from. Copy chosen DOIs into the sources_doi column of data/assays.csv (separate several with ;).
--count-only prints the number of matching papers: a rough input for n_literature_using, NOT a precise count.
"""
import argparse, csv, json, os, urllib.parse, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def search(query, size):
    url = API + "?" + urllib.parse.urlencode({"query": query, "format": "json", "resultType": "lite", "pageSize": size, "sort": "CITED desc"})
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "assay-builder-helper"}), timeout=60) as r:
        return json.load(r)


def parse(payload):
    return [{"doi": x.get("doi", ""), "pmid": x.get("pmid", ""), "title": x.get("title", ""), "year": x.get("pubYear", ""),
             "journal": x.get("journalTitle", ""), "cited_by": x.get("citedByCount", ""), "open_access": x.get("isOpenAccess", "")}
            for x in payload.get("resultList", {}).get("result", [])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("query"); ap.add_argument("--max", type=int, default=25); ap.add_argument("--count-only", action="store_true")
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "literature_candidates.csv"))
    a = ap.parse_args()
    data = search(a.query, 1 if a.count_only else a.max)
    print("Europe PMC reports %s matching papers." % data.get("hitCount"))
    if a.count_only:
        return
    rows = parse(data)
    with open(a.out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["doi", "pmid", "title", "year", "journal", "cited_by", "open_access"]); w.writeheader(); w.writerows(rows)
    print("Wrote %d rows to %s" % (len(rows), a.out))


if __name__ == "__main__":
    main()
