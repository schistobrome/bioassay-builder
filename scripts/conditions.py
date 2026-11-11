#!/usr/bin/env python3
"""
Assay conditions (pH, temperature, buffer, detergent, DMSO, ATP, incubation time ...) stated in the assay texts.

    python scripts/conditions.py                 # all summarized proteins (a few minutes; never uses the internet)
    python scripts/conditions.py --redo          # recompute everything
    python scripts/build_database.py --aggregate-only --merge --replace-m2     # then summarize them per assay type

Reads the same texts as plate_formats.py. For every campaign it records the FIRST value the text states for each
parameter (a ChEMBL copy of a PubChem deposit counts with the PubChem original) in build_cache/conditions/<accession>.json.
build_database.py then writes median and spread per assay type to data/condition_stats.csv.
Not covered: BindingDB (its pH and temperature columns are not read yet). Texts often mention several buffers
(reaction, detection, wash); the first stated value is used, so treat the numbers as typical, not exact.
"""
import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_database as bd      # noqa: E402
from build_database import m1, m2   # noqa: E402


def campaigns_with_conditions(key):
    records, _ = m2.build_records(key, m1.DATA, None)
    if not records:
        return {}
    original = {}
    for r in records:
        if (r.get("source") or "").startswith("PubChem/") and r.get("pubchem_aid"):
            original[str(r["pubchem_aid"]).strip()] = r.get("campaign") or r["id"]
    out = {}
    for r in records:
        camp = r.get("campaign") or r["id"]
        aid = str(r.get("pubchem_aid") or "").strip()
        if r.get("source") == "ChEMBL" and aid in original and original[aid] != camp:
            camp = original[aid]
        found = bd.conditions_in(r.get("text"))
        slot = out.setdefault(camp, {})
        for k, v in found.items():
            slot.setdefault(k, v)                     # first record that states it wins
    return {k: v for k, v in out.items() if v}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--redo", action="store_true", help="recompute proteins that already have a conditions file")
    ap.add_argument("--max-minutes", type=float, default=0, help="stop cleanly after this long (0 = no limit)")
    args = ap.parse_args(argv)
    m1.DATA = os.path.join(bd.CACHE, "m1_data")
    os.makedirs(bd.CONDS_DIR, exist_ok=True)
    files = sorted(f for f in os.listdir(bd.SUMMARY) if f.endswith(".json")) if os.path.isdir(bd.SUMMARY) else []
    started, done, skipped, stated = time.time(), 0, 0, {}
    for n, f in enumerate(files, 1):
        s = bd.read_json(os.path.join(bd.SUMMARY, f))
        if not s:
            continue
        path = os.path.join(bd.CONDS_DIR, f)
        if os.path.exists(path) and not args.redo:
            skipped += 1
            continue
        if args.max_minutes and (time.time() - started) / 60 > args.max_minutes:
            print("Stopped at the time limit; run again to continue.")
            break
        camps = campaigns_with_conditions(s["gene"])
        for k in camps.values():
            for p in k:
                stated[p] = stated.get(p, 0) + 1
        bd.write_json(path, {"acc": s["acc"], "gene": s["gene"], "built": bd.now(), "campaigns": camps})
        done += 1
        if done % 250 == 0:
            print("  %d proteins done (%d of %d)" % (done, n, len(files)))
    print("Conditions: %d proteins done, %d already done." % (done, skipped))
    print("Campaigns stating each parameter (this run): " + ", ".join("%s %d" % (bd.COND_LABEL[p], stated.get(p, 0))
                                                                    for p in bd.COND_NUMERIC + bd.COND_CATEGORICAL))
    print("Now run: python scripts/build_database.py --aggregate-only --merge --replace-m2")


if __name__ == "__main__":
    main()
