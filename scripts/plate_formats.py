#!/usr/bin/env python3
"""
Plate formats (96 / 384 / 1536 / low-density) stated in the assay texts of every protein build_database.py summarized.

    python scripts/plate_formats.py                 # all summarized proteins (a few minutes; never uses the internet)
    python scripts/plate_formats.py --redo          # recompute everything
    python scripts/build_database.py --aggregate-only --merge --replace-m2     # then fold them into the usage rows

It reads the assay texts already on disk (ChEMBL descriptions and the PubChem descriptions and protocols),
looks for statements such as "384-well" or "1536 well plates", ignores compound / source / stock plates, and writes
build_cache/plates/<accession>.json with the formats per campaign (a ChEMBL copy of a PubChem deposit counts with the
PubChem original, as in the usage statistics). Not covered: BindingDB assays (their text rarely names a plate).
"""
import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_database as bd      # noqa: E402
from build_database import m1, m2   # noqa: E402


def campaigns_with_plates(key):
    records, _ = m2.build_records(key, m1.DATA, None)
    if not records:
        return {}, 0
    original = {}
    for r in records:
        if (r.get("source") or "").startswith("PubChem/") and r.get("pubchem_aid"):
            original[str(r["pubchem_aid"]).strip()] = r.get("campaign") or r["id"]
    out, with_text = {}, 0
    for r in records:
        camp = r.get("campaign") or r["id"]
        aid = str(r.get("pubchem_aid") or "").strip()
        if r.get("source") == "ChEMBL" and aid in original and original[aid] != camp:
            camp = original[aid]
        found = bd.plates_in(r.get("text"))
        out.setdefault(camp, set()).update(found)
        with_text += bool(found)
    return {k: sorted(v) for k, v in out.items()}, with_text


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--redo", action="store_true", help="recompute proteins that already have a plates file")
    ap.add_argument("--max-minutes", type=float, default=0, help="stop cleanly after this long (0 = no limit)")
    args = ap.parse_args(argv)
    m1.DATA = os.path.join(bd.CACHE, "m1_data")          # where build_database.py keeps the downloaded assay texts
    os.makedirs(bd.PLATES_DIR, exist_ok=True)
    files = sorted(f for f in os.listdir(bd.SUMMARY) if f.endswith(".json")) if os.path.isdir(bd.SUMMARY) else []
    started, done, skipped, empty = time.time(), 0, 0, 0
    totals = {"96": 0, "384": 0, "1536": 0, "low": 0}
    for n, f in enumerate(files, 1):
        s = bd.read_json(os.path.join(bd.SUMMARY, f))
        if not s:
            continue
        path = os.path.join(bd.PLATES_DIR, f)
        if os.path.exists(path) and not args.redo:
            skipped += 1
            continue
        if args.max_minutes and (time.time() - started) / 60 > args.max_minutes:
            print("Stopped at the time limit; run again to continue.")
            break
        camps, with_text = campaigns_with_plates(s["gene"])
        for fmts in camps.values():
            for p in fmts:
                totals[p] = totals.get(p, 0) + 1
        empty += not any(camps.values())
        bd.write_json(path, {"acc": s["acc"], "gene": s["gene"], "built": bd.now(), "campaigns": camps})
        done += 1
        if done % 250 == 0:
            print("  %d proteins done (%d of %d)" % (done, n, len(files)))
    print("Plate formats: %d proteins done, %d already done, %d with no plate format stated anywhere." % (done, skipped, empty))
    print("Campaigns stating each format (this run): 96-well %(96)d, 384-well %(384)d, 1536-well %(1536)d, 6-48-well %(low)d" % totals)
    print("Now run: python scripts/build_database.py --aggregate-only --merge --replace-m2")


if __name__ == "__main__":
    main()
