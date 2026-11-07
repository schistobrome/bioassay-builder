"""
M1b: download the original PubChem descriptions that ChEMBL assays point to.

ChEMBL sometimes holds only a one-line title for an assay that was deposited in
PubChem. M2 writes the list of those PubChem numbers (AIDs) into
m2_data/<GENE>_original_pubchem_aids.txt. This script downloads each one's full
description and protocol into m1_data/pubchem_raw (the same place M1 uses), and
skips any it already has.

How to run, after M2 has run once:
    .\\.venv\\Scripts\\python.exe m1b_fetch_original_aids.py
    (or with your own numbers:  m1b_fetch_original_aids.py 651768 720497)

Then run m2_label.py again. It reads these files and merges the text in.
It only READS public data. Uses only Python's built-in libraries.
"""

import glob
import json
import os
import sys
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "m1_data", "pubchem_raw")
HEADERS = {"User-Agent": "assay-planner-m1b/0.1 (personal side project)", "Accept": "application/json"}


def fetch(aid):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/description/JSON" % aid
    last = None
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504):
                raise
            last = error
        except (urllib.error.URLError, TimeoutError) as error:
            last = error
        time.sleep(3 * (attempt + 1))
    raise last


def main():
    wanted = [a for a in sys.argv[1:] if a.isdigit()]
    if not wanted:
        for path in glob.glob(os.path.join(HERE, "m2_data", "*_original_pubchem_aids.txt")):
            with open(path) as handle:
                wanted += [line.strip() for line in handle if line.strip().isdigit()]
    wanted = sorted(set(wanted), key=int)
    if not wanted:
        print("Nothing to download. Run m2_label.py first, or give AID numbers.")
        return
    os.makedirs(RAW, exist_ok=True)
    done = failed = 0
    for aid in wanted:
        path = os.path.join(RAW, "AID_%s.json" % aid)
        if os.path.exists(path):
            print("  AID %s: already have it" % aid)
            continue
        try:
            data = fetch(aid)
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(data, handle)
            name = data["PC_AssayContainer"][0]["assay"]["descr"].get("name", "")
            print("  AID %s: saved (%s)" % (aid, " ".join(str(name).split())[:70]))
            done += 1
        except Exception as error:
            print("  AID %s: FAILED (%s: %s)" % (aid, type(error).__name__, error))
            failed += 1
        time.sleep(0.4)
    print("Done: %d downloaded, %d failed. Now run m2_label.py again." % (done, failed))


if __name__ == "__main__":
    main()
