#!/usr/bin/env python3
"""
Turn the downloaded PubChem BioAssay zips into the files build_database.py reads (no internet needed).

    python scripts/pubchem_local.py                       # after pubchem_plan.py and the zip download
    python scripts/build_database.py --pubchem-local ...  # the build then makes NO PubChem calls

For every assay listed in build_cache/pubchem_bulk/needed_aids.csv it opens the matching zip in
data/public/pubchem/JSON/, reads <aid>.json.gz, and writes
  build_cache/m1_data/pubchem_raw/AID_<aid>.json     the description in the format the labeler reads
  build_cache/pubchem_bulk/csv/<protein key>.csv      one table per protein (same columns the API fetch made)
Resumable and quick to rerun (existing files are kept unless --redo). Standard library only.
"""
import argparse
import collections
import csv
import gzip
import json
import os
import re
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CACHE = os.path.join(ROOT, "build_cache")
BULK = os.path.join(CACHE, "pubchem_bulk")
RAW = os.path.join(CACHE, "m1_data", "pubchem_raw")
COLUMNS = ["aid", "source", "title", "found_by", "summary_text", "full_description", "protocol", "targets"]


def safe(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def text_of(value):
    if isinstance(value, list):
        return "\n".join(str(v) for v in value if v is not None)
    return "" if value is None else str(value)


def clip(text, n=1500):
    text = re.sub(r"\s+", " ", text or "").strip()
    return text[:n]


def zip_name(aid, width=1000):
    lo = (aid - 1) // width * width + 1
    return "%07d_%07d" % (lo, lo + width - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-dir", default="data/public/pubchem/JSON")
    ap.add_argument("--redo", action="store_true")
    args = ap.parse_args()
    need_path = os.path.join(BULK, "needed_aids.csv")
    if not os.path.exists(need_path):
        sys.exit("Run pubchem_plan.py first (build_cache/pubchem_bulk/needed_aids.csv is missing).")
    with open(need_path, newline="", encoding="utf-8-sig") as h:
        needed = list(csv.DictReader(h))
    with open(os.path.join(CACHE, "uniprot_human_reviewed.csv"), newline="", encoding="utf-8-sig") as h:
        proteins = {r["accession"]: r for r in csv.DictReader(h)}
    counts = collections.Counter(p["gene"] for p in proteins.values() if p["gene"])

    def key_of(acc):
        p = proteins[acc]
        return safe(p["gene"]) if p["gene"] and counts[p["gene"]] == 1 else safe((p["gene"] or "NA") + "_" + acc)

    by_zip = collections.defaultdict(set)
    for r in needed:
        by_zip[zip_name(int(r["aid"]))].add(int(r["aid"]))
    os.makedirs(RAW, exist_ok=True)
    info, missing_zip, missing_file, bad = {}, 0, 0, 0
    for z_i, (zname, aids) in enumerate(sorted(by_zip.items()), 1):
        path = os.path.join(args.json_dir, zname + ".zip")
        if not os.path.exists(path):
            missing_zip += len(aids)
            print("  [%d/%d] %s not downloaded (%d assays skipped)" % (z_i, len(by_zip), zname, len(aids)))
            continue
        with zipfile.ZipFile(path) as z:
            names = set(z.namelist())
            for aid in sorted(aids):
                member = "%s/%d.json.gz" % (zname, aid)
                if member not in names:
                    missing_file += 1
                    continue
                try:
                    descr = json.loads(gzip.decompress(z.read(member)))["PC_AssaySubmit"]["assay"]["descr"]
                except (KeyError, ValueError, OSError, TypeError):
                    bad += 1
                    continue
                target = os.path.join(RAW, "AID_%d.json" % aid)
                if args.redo or not os.path.exists(target) or os.path.getsize(target) == 0:
                    tmp = target + ".tmp"
                    with open(tmp, "w", encoding="utf-8") as f:
                        json.dump({"PC_AssayContainer": [{"assay": {"descr": descr}}]}, f)
                    os.replace(tmp, target)
                targets = "; ".join(t["name"] for t in (descr.get("target") or []) if isinstance(t, dict) and t.get("name"))
                src = ((descr.get("aid_source") or {}).get("db") or {}).get("name", "")
                info[aid] = {"title": text_of(descr.get("name")), "source": src,
                             "summary": clip(text_of(descr.get("description")), 1500),
                             "full": clip(text_of(descr.get("description")), 1500),
                             "protocol": clip(text_of(descr.get("protocol")), 1500), "targets": targets}
        print("  [%d/%d] %s: %d assays read" % (z_i, len(by_zip), zname, len(aids)))

    per_key = collections.defaultdict(list)
    for r in needed:
        aid = int(r["aid"])
        if aid not in info or r["accession"] not in proteins:
            continue
        i = info[aid]
        per_key[key_of(r["accession"])].append(
            {"aid": aid, "source": r["source"] or i["source"], "title": i["title"] or r["name"],
             "found_by": "local bulk", "summary_text": i["summary"], "full_description": i["full"],
             "protocol": i["protocol"], "targets": i["targets"]})
    out = os.path.join(BULK, "csv")
    os.makedirs(out, exist_ok=True)
    for key, rows in per_key.items():
        with open(os.path.join(out, key + ".csv"), "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=COLUMNS)
            w.writeheader()
            w.writerows(rows)
    print("Wrote %d protein tables (%d assay rows). Not downloaded: %d, not found in zips: %d, unreadable: %d."
          % (len(per_key), sum(len(v) for v in per_key.values()), missing_zip, missing_file, bad))
    print("Next: python scripts/build_database.py --pubchem-local --max-minutes 480 --pause 1.0 --timeout 90")


if __name__ == "__main__":
    main()
