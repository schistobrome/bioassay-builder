"""
M2b: download the title and abstract of every paper behind the ChEMBL assays.

One ChEMBL call fetches up to 50 papers, and every paper is saved once as
m1_data/chembl_docs/<ID>.json, so later runs call nothing for papers you have.

How to run (after M1):
    .\\.venv\\Scripts\\python.exe m2b_fetch_abstracts.py
    (or only one gene:  m2b_fetch_abstracts.py egfr)
Then run m2_label.py again. It reads these files automatically.
It only READS public data. Uses only Python's built-in libraries.
"""

import csv
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
M1 = os.path.join(HERE, "m1_data")
OUT = os.path.join(M1, "chembl_docs")
HEADERS = {"User-Agent": "assay-planner-m2b/0.1 (personal side project)", "Accept": "application/json"}
BASE = "https://www.ebi.ac.uk/chembl/api/data/document.json"
KEEP = ("document_chembl_id", "title", "abstract", "pubmed_id", "doi", "journal", "year", "doc_type")
BATCH = 50


def http_get_text(url, timeout=90):
    request = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def fetch_batch(ids):
    url = BASE + "?" + urllib.parse.urlencode({"document_chembl_id__in": ",".join(ids), "limit": 100})
    last = None
    for attempt in range(3):
        try:
            return json.loads(http_get_text(url)).get("documents", [])
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504):
                raise
            last = error
        except (urllib.error.URLError, TimeoutError) as error:
            last = error
        time.sleep(3 * (attempt + 1))
    raise last


def wanted_ids(genes):
    ids = set()
    for gene in genes:
        path = os.path.join(M1, "chembl_%s_assays.csv" % gene)
        if not os.path.exists(path):
            print("  No %s - run M1 first." % path)
            continue
        with open(path, newline="", encoding="utf-8-sig") as handle:
            ids.update(r["document"] for r in csv.DictReader(handle) if r.get("document"))
    return sorted(ids)


def main():
    genes = [g.upper() for g in sys.argv[1:]] or ["WRN", "EGFR"]
    ids = wanted_ids(genes)
    os.makedirs(OUT, exist_ok=True)
    todo = [i for i in ids if not os.path.exists(os.path.join(OUT, i + ".json"))]
    print("Papers needed: %d. Already saved: %d. To download: %d." % (len(ids), len(ids) - len(todo), len(todo)))
    saved = failed = 0
    for start in range(0, len(todo), BATCH):
        part = todo[start:start + BATCH]
        try:
            docs = fetch_batch(part)
        except Exception as error:
            print("  batch %d FAILED (%s: %s). Run again to retry." % (start // BATCH + 1, type(error).__name__, error))
            failed += len(part)
            continue
        got = set()
        for doc in docs:
            record = {k: doc.get(k) for k in KEEP}
            with open(os.path.join(OUT, "%s.json" % doc["document_chembl_id"]), "w", encoding="utf-8") as handle:
                json.dump(record, handle)
            got.add(doc["document_chembl_id"])
            saved += 1
        for missing in set(part) - got:     # remember papers ChEMBL has nothing for
            with open(os.path.join(OUT, missing + ".json"), "w", encoding="utf-8") as handle:
                json.dump({"document_chembl_id": missing}, handle)
        print("  %d of %d papers done" % (min(start + BATCH, len(todo)), len(todo)))
        time.sleep(0.4)
    print("Done: %d saved, %d failed. Now run m2_label.py again." % (saved, failed))


if __name__ == "__main__":
    main()
