#!/usr/bin/env python3
"""
Decide which PubChem BioAssay description zips you really need (no internet).

    python scripts/pubchem_plan.py --extras data/public/pubchem/Extras

Reads bioassays.tsv.gz (+ RevokedAids.gz) and keeps assays that
  * have a UniProt target that is in build_cache/uniprot_human_reviewed.csv (your human proteins),
  * are not revoked, not RNAi (nucleotide) screens,
  * are not copies of data you already have: depositor ChEMBL (comes from the ChEMBL download)
    and BindingDB (comes from the BindingDB download).
Writes, in build_cache/pubchem_bulk/:
  needed_aids.csv        one row per assay and target protein
  needed_json_zips.txt   names of the JSON/ zip files that contain those assays
and prints the exact download command. Standard library only.
"""
import argparse
import collections
import csv
import gzip
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "build_cache", "pubchem_bulk")
ACC = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")
COPY_SOURCES = {"chembl", "bindingdb"}


def zip_name(aid, width=1000):
    lo = (aid - 1) // width * width + 1
    return "%07d_%07d.zip" % (lo, lo + width - 1)


def pick(header, *needles):
    for i, h in enumerate(header):
        if all(n in h.lower() for n in needles):
            return i
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extras", default="data/public/pubchem/Extras")
    ap.add_argument("--reviewed", default=os.path.join(ROOT, "build_cache", "uniprot_human_reviewed.csv"))
    ap.add_argument("--keep-copies", action="store_true", help="also keep ChEMBL / BindingDB depositor assays")
    args = ap.parse_args()
    path = os.path.join(args.extras, "bioassays.tsv.gz")
    if not os.path.exists(path):
        sys.exit("Not found: %s (download the Extras files first)" % path)
    reviewed = set()
    if os.path.exists(args.reviewed):
        with open(args.reviewed, newline="", encoding="utf-8-sig") as h:
            reviewed = {r["accession"] for r in csv.DictReader(h)}
    if not reviewed:
        sys.exit("Missing %s: run build_database.py once so the human protein list exists." % args.reviewed)
    revoked = set()
    rpath = os.path.join(args.extras, "RevokedAids.gz")
    if os.path.exists(rpath):
        with gzip.open(rpath, "rt", encoding="utf-8", errors="replace") as h:
            for line in h:
                m = re.match(r"\s*(\d+)", line)
                if m:
                    revoked.add(int(m.group(1)))
    print("Human proteins: %d | revoked assays: %d" % (len(reviewed), len(revoked)))

    rows, seen, why = [], 0, collections.Counter()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as h:
        reader = csv.reader(h, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader)
        c = {"aid": pick(header, "aid"), "name": pick(header, "bioassay name"), "src": pick(header, "source name"),
             "sub": pick(header, "substance type"), "out": pick(header, "outcome type"),
             "types": pick(header, "bioassay types"), "uni": pick(header, "uniprot")}
        if None in (c["aid"], c["src"], c["uni"]):
            print("HEADER SEEN:", header)
            sys.exit("Could not find the AID / Source Name / UniProt columns; send me this header.")
        for row in reader:
            seen += 1
            try:
                aid = int(row[c["aid"]])
            except (ValueError, IndexError):
                continue
            accs = {m for m in ACC.findall(row[c["uni"]] if c["uni"] < len(row) else "")} & reviewed
            if not accs:
                why["no human reviewed target"] += 1
                continue
            if aid in revoked:
                why["revoked"] += 1
                continue
            src = row[c["src"]].strip()
            if not args.keep_copies and src.lower() in COPY_SOURCES:
                why["copy (depositor %s)" % src] += 1
                continue
            if c["sub"] is not None and "nucleotide" in row[c["sub"]].lower() and "small" not in row[c["sub"]].lower():
                why["RNAi / nucleotide screen"] += 1
                continue
            for a in sorted(accs):
                rows.append({"aid": aid, "accession": a, "source": src,
                             "outcome": row[c["out"]] if c["out"] is not None else "",
                             "types": row[c["types"]] if c["types"] is not None else "",
                             "name": row[c["name"]][:150] if c["name"] is not None else ""})
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "needed_aids.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["aid", "accession", "source", "outcome", "types", "name"])
        w.writeheader()
        w.writerows(rows)
    aids = sorted({r["aid"] for r in rows})
    zips = sorted({zip_name(a) for a in aids})
    with open(os.path.join(OUT, "needed_json_zips.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(zips) + "\n")
    per_protein = collections.Counter(r["accession"] for r in rows)
    print("Read %d assays in bioassays.tsv.gz" % seen)
    print("Left out: " + "; ".join("%d %s" % (v, k) for k, v in why.most_common()))
    print("KEEP: %d assays on %d human proteins -> %d description zips needed" % (len(aids), len(per_protein), len(zips)))
    print("Busiest proteins: " + ", ".join("%s %d" % kv for kv in per_protein.most_common(8)))
    print("Top depositors kept: " + "; ".join("%s %d" % kv for kv in collections.Counter(r["source"] for r in rows).most_common(8)))
    print("\nWrote build_cache/pubchem_bulk/needed_aids.csv and needed_json_zips.txt. Next:\n"
          "  python scripts/download_bulk.py --url https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/JSON/ "
          "--dest data/public/pubchem/JSON --only-from build_cache/pubchem_bulk/needed_json_zips.txt")


if __name__ == "__main__":
    main()
