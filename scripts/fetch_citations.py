#!/usr/bin/env python3
"""Fill data/citations.csv (title, year, journal, authors) for every DOI listed in data/assays.csv sources_doi.

    python scripts/fetch_citations.py
Uses the free Crossref API. DOIs already in citations.csv are skipped.
"""
import csv, json, os, re, urllib.parse, urllib.request
from datetime import date

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIELDS = ["doi", "pmid", "title", "year", "journal", "authors", "type", "url", "fetched_on"]


def crossref(doi):
    req = urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(doi), headers={"User-Agent": "assay-builder-helper"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["message"]


def to_row(doi, m):
    year = ((m.get("issued") or {}).get("date-parts") or [[""]])[0][0]
    authors = "; ".join("%s %s" % (a.get("family", ""), (a.get("given", "") or "")[:1]) for a in m.get("author", [])[:6])
    return {"doi": doi, "pmid": "", "title": (m.get("title") or [""])[0], "year": year, "journal": (m.get("container-title") or [""])[0],
            "authors": authors, "type": "", "url": "https://doi.org/" + doi, "fetched_on": date.today().isoformat()}


def main():
    cpath = os.path.join(ROOT, "data", "citations.csv")
    have = {r["doi"].lower() for r in csv.DictReader(open(cpath, encoding="utf-8-sig"))}
    dois = []
    for r in csv.DictReader(open(os.path.join(ROOT, "data", "assays.csv"), encoding="utf-8-sig")):
        for d in re.split(r"[;|]", r.get("sources_doi") or ""):
            d = re.sub(r"^https?://(dx\.)?doi\.org/", "", d.strip())
            if d and d.lower() not in have and d not in dois:
                dois.append(d)
    new = []
    for d in dois:
        try:
            new.append(to_row(d, crossref(d))); print("ok  ", d)
        except Exception as e:
            print("FAIL", d, e)
    if new:
        with open(cpath, "a", newline="", encoding="utf-8-sig") as f:
            csv.DictWriter(f, fieldnames=FIELDS).writerows(new)
    print("Added %d citation(s)." % len(new))


if __name__ == "__main__":
    main()
