#!/usr/bin/env python3
"""
BindingDB -> labeled campaigns for the usage statistics.

    python scripts/bindingdb_label.py --folder data/public/bindingdb
    python scripts/bindingdb_label.py --folder data/public/bindingdb --rules-only --max-minutes 60
    python scripts/build_database.py --aggregate-only --merge          # then fold the new campaigns into the usage rows

What it does (it never touches the internet):
  1. One pass over BindingDB_All*.tsv (about 4 minutes the first time, then cached in build_cache/bdb/agg.json):
     groups the measurements per protein and per assay, remembers the data source, patent number, PMID,
     affinity types and whether the target name names a mutant.
  2. For every protein that build_database.py already summarized (or --all): takes the assay descriptions,
     labels them with the same m2 labeler (function, technology, cell-free / cell-based, construct range),
     and writes build_cache/summary_bdb/<accession>.json (campaigns in the same format as ChEMBL / PubChem).
Direction (inhibition / activation) comes from the start of the assay description plus the affinity types,
as in bindingdb_local.py, and replaces the labeler's whole-text guess.

NOT double counted: rows BindingDB copied from ChEMBL or PubChem are skipped. Patent rows get the campaign
"patent:<number>", papers "pmid:<id>" (or "bdb:<entry>"). BindingDB's own literature curation can overlap
papers ChEMBL also has; use --no-literature to count patents and other sources only.
Resumable: proteins already done are skipped; --max-minutes stops cleanly. Standard library plus build_database.py.
"""
import argparse
import collections
import contextlib
import csv
import html
import io
import os
import re
import sys
import time
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_database as bd          # noqa: E402
import bindingdb_local as bl         # noqa: E402

m1, m2 = bd.m1, bd.m2
CACHE = bd.CACHE
BDB = os.path.join(CACHE, "bdb")
SUMMARY_BDB = os.path.join(CACHE, "summary_bdb")
LABELS_BDB = os.path.join(CACHE, "m2_data_bdb")
COPY_SOURCES = ("chembl", "pubchem")
EFFECT_MAP = {"inhibition": "inhibition", "inhibition_ic50_only": "inhibition", "activation": "activation",
              "activation_ec50_only": "activation", "mixed": "mixed", "binding_only": "unspecified",
              "unclear": "unspecified"}
MUT = re.compile(r"\[[^\]]*\b[A-Z]\d{1,4}[A-Z*]\b[^\]]*\]")


def aggregate_table(folder, reviewed, cache_path):
    main_path = bl.find_file(folder, "BindingDB_All*.tsv") or bl.find_file(folder, "BindingDB*.tsv",
                                                                          exclude=("assays", "rsid", "uniprot"))
    rsid_path = bl.find_file(folder, "*rsid*eaid*") or bl.find_file(folder, "*rsid*")
    if not main_path or not rsid_path:
        sys.exit("Need BindingDB_All*.tsv and the rsid_eaids file inside %s" % folder)
    meta = {"size": os.path.getsize(main_path), "mtime": int(os.path.getmtime(main_path)), "reviewed": len(reviewed)}
    cached = bd.read_json(cache_path)
    if cached and cached.get("meta") == meta:
        print("Using cached BindingDB table (%s)" % cache_path)
        return cached["agg"]
    rs_to_assay = {}
    h, r = bl.read_tsv(rsid_path)
    next(r, None)
    for row in r:
        if len(row) >= 2:
            rs_to_assay[row[0].strip()] = row[1].strip()
    h.close()
    handle, reader = bl.read_tsv(main_path)
    header = next(reader)
    c_acc = bl.pick(header, "uniprot", "primary")
    c_rs = bl.pick(header, "reactant_set") or bl.pick(header, "reaction", "set")
    c_name = bl.pick(header, "target name")
    c_src = bl.pick(header, "curation") or bl.pick(header, "datasource")
    c_pat = bl.pick(header, "patent")
    c_pmid = bl.pick(header, "pmid")
    c_aff = {t: bl.pick(header, t.lower() + " (nm)") for t in bl.AFFINITY}
    if not c_acc or not c_rs:
        print("HEADER SEEN:")
        for i, x in enumerate(header):
            print("  %3d %s" % (i, x))
        sys.exit("Could not find the UniProt or reaction-set columns; send me this header list.")
    print("Reading %s (patent col %s, pmid col %s)" % (main_path, bool(c_pat), bool(c_pmid)))
    agg = collections.defaultdict(dict)
    width, t0, rows = len(header), time.time(), 0
    for row in reader:
        rows += 1
        if rows % 500000 == 0:
            print("  %d rows (%.0fs)" % (rows, time.time() - t0))
        if len(row) < width:
            row += [""] * (width - len(row))
        accs = set()
        for i in c_acc:
            accs.update(m.upper() for m in bl.ACC.findall(row[i]))
        keep = accs & reviewed
        if not keep:
            continue
        rs = row[c_rs[0]].strip()
        key = rs_to_assay.get(rs) or "rs" + rs
        types = {t for t, cols in c_aff.items() if any(row[i].strip() for i in cols)}
        src = row[c_src[0]].strip() if c_src else ""
        pat = row[c_pat[0]].strip() if c_pat else ""
        pmid = row[c_pmid[0]].strip() if c_pmid else ""
        mut = bool(MUT.search(row[c_name[0]])) if c_name else False
        for a in keep:
            e = agg[a].get(key)
            if e is None:
                agg[a][key] = [1, "".join(sorted(t[0] if t != "Kd" else "D" for t in types)), src, pat, pmid, int(mut), ",".join(sorted(types))]
            else:
                e[0] += 1
                e[5] |= int(mut)
                merged = set(e[6].split(",")) | types
                e[6] = ",".join(sorted(x for x in merged if x))
    handle.close()
    print("Read %d rows; %d proteins, %d assays (%.0fs)" % (rows, len(agg), sum(len(v) for v in agg.values()), time.time() - t0))
    bd.write_json(cache_path, {"meta": meta, "agg": agg})
    return agg


def campaign_of(key, pat, pmid):
    if pat:
        return "patent:" + re.sub(r"[^A-Za-z0-9]", "", pat).upper()
    if pmid:
        return "pmid:" + pmid
    return "bdb:" + key.split("_")[0]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", required=True)
    ap.add_argument("--all", action="store_true", help="every reviewed protein in BindingDB, not only those build_database summarized")
    ap.add_argument("--genes", help="comma-separated gene symbols or accessions")
    ap.add_argument("--rules-only", action="store_true", help="synonym rules only (no meaning model)")
    ap.add_argument("--no-literature", action="store_true", help="skip BindingDB's own literature curation (may overlap ChEMBL)")
    ap.add_argument("--max-assays-per-protein", type=int, default=3000)
    ap.add_argument("--max-minutes", type=float)
    ap.add_argument("--redo", action="store_true", help="label everything again")
    args = ap.parse_args(argv)

    proteins = bd.read_csv(os.path.join(CACHE, "uniprot_human_reviewed.csv"))
    if not proteins:
        sys.exit("build_cache/uniprot_human_reviewed.csv is missing: run build_database.py once first.")
    reviewed = {p["accession"] for p in proteins}
    counts = collections.Counter(p["gene"] for p in proteins if p["gene"])
    info = {p["accession"]: p for p in proteins}
    agg = aggregate_table(args.folder, reviewed, os.path.join(BDB, "agg.json"))
    assay_text = {}
    path = bl.find_file(args.folder, "BindingDB_Assays*")
    if not path:
        sys.exit("BindingDB_Assays*.tsv not found")
    h, r = bl.read_tsv(path)
    next(r, None)
    for row in r:
        if len(row) >= 3:
            assay_text["%s_%s" % (row[0].strip(), row[1].strip())] = " ".join(x.strip() for x in row[2:] if x.strip())
    h.close()

    done_main = {f[:-5] for f in os.listdir(bd.SUMMARY)} if os.path.isdir(bd.SUMMARY) else set()
    if args.genes:
        names = {g.strip().upper() for g in args.genes.split(",") if g.strip()}
        targets = [a for a in agg if a in names or info[a]["gene"].upper() in names]
    elif args.all:
        targets = list(agg)
    else:
        targets = [a for a in agg if bd.safe(a) in done_main]
    targets.sort(key=lambda a: -len(agg[a]))
    print("%d proteins to label (%s)" % (len(targets), "all" if args.all else "those already summarized by build_database.py"))

    added = bd.ensure_terms()
    terms = m2.load_terms()
    term_hash = bd.file_hash(m2.TERM_TABLE) + ("r" if args.rules_only else "m") + \
        ("L" if args.no_literature else "") + str(args.max_assays_per_protein)
    embedder = bd.EmbedderBox(not args.rules_only)
    bd.log.handle = None
    m2.say = lambda *a, **k: None
    for folder in (SUMMARY_BDB, LABELS_BDB, os.path.join(BDB, "state")):
        os.makedirs(folder, exist_ok=True)
    deadline = time.time() + args.max_minutes * 60 if args.max_minutes else None
    n_done = n_skip = 0
    try:
        for n, acc in enumerate(targets, 1):
            state_path = os.path.join(BDB, "state", bd.safe(acc) + ".json")
            st = bd.read_json(state_path, {})
            if not args.redo and st.get("term_hash") == term_hash and st.get("built_from") == len(agg[acc]):
                n_skip += 1
                continue
            if deadline and time.time() >= deadline:
                print("  time limit reached; rerun the same command to continue")
                break
            row = info[acc]
            key = bd.safe(row["gene"]) if row["gene"] and counts[row["gene"]] == 1 else bd.safe((row["gene"] or "NA") + "_" + acc)
            items = sorted(agg[acc].items(), key=lambda kv: -kv[1][0])
            records, skipped = [], collections.Counter()
            for akey, (cnt, _m, src, pat, pmid, mut, types) in items:
                low = src.lower()
                if low.startswith(COPY_SOURCES):
                    skipped["copy of ChEMBL/PubChem"] += 1
                    continue
                if args.no_literature and "literature" in low:
                    skipped["literature (--no-literature)"] += 1
                    continue
                text = html.unescape(assay_text.get(akey, "")).strip()
                if not text:
                    skipped["no description"] += 1
                    continue
                if len(records) >= args.max_assays_per_protein:
                    skipped["over per-protein cap"] += 1
                    continue
                records.append({"id": "bdb:" + akey, "source": "BindingDB/" + (src or "?"),
                                "campaign": campaign_of(akey, pat, pmid), "text_level": "BindingDB assay description",
                                "text": text[:3000], "title": text[:150], "chembl_type": "", "bao": "", "cell_type": "",
                                "mutation": "yes" if mut else "", "pubchem_aid": "", "abstract": "",
                                "_types": set(t for t in types.split(",") if t)})
            if records:
                length = m2.protein_length(os.path.join(CACHE, "m1_data"), acc)
                implied = m2.implied_function(terms, os.path.join(CACHE, "m1_data"), acc)
                emb = embedder.get()
                for i in range(0, len(records), 300):
                    if deadline and time.time() >= deadline + 600:        # allow a protein to finish, but not forever
                        break
                    chunk = records[i:i + 300]
                    with contextlib.redirect_stdout(io.StringIO()):
                        m2.label_records(chunk, terms, emb, m2.construct_patterns(key), length, implied)
                for rec in records:
                    rec["effect"] = EFFECT_MAP[bl.effect_of(rec["text"], rec["_types"])]
            bd.write_csv(os.path.join(LABELS_BDB, key + "_bdb_labeled.csv"), records, bd.LABEL_COLUMNS)
            campaigns, merged = bd.dedupe_campaigns(records)
            bd.write_json(os.path.join(SUMMARY_BDB, bd.safe(acc) + ".json"),
                          {"acc": acc, "n_assays": len(records), "campaigns": campaigns, "skipped": dict(skipped),
                           "built": bd.now()})
            bd.write_json(state_path, {"term_hash": term_hash, "built_from": len(agg[acc])})
            n_done += 1
            print("  [%d/%d] %-10s %6d assays -> %5d campaigns  (skipped: %s)"
                  % (n, len(targets), row["gene"] or acc, len(records), len(campaigns),
                     ", ".join("%d %s" % (v, k) for k, v in skipped.items()) or "none"))
            if n_done % 25 == 0:
                embedder.save()
    except KeyboardInterrupt:
        print("\nStopped by you. Rerun to continue.")
    finally:
        embedder.save()
    print("Done: %d proteins labeled, %d already done. Now run:  python scripts/build_database.py --aggregate-only --merge"
          % (n_done, n_skip))


if __name__ == "__main__":
    main()
