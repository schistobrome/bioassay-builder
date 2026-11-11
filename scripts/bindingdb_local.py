#!/usr/bin/env python3
"""
Summarize a local BindingDB download per UniProt protein (no internet needed).

    python scripts/bindingdb_local.py --folder C:\\path\\to\\bindingdb_unzipped
    python scripts/bindingdb_local.py --folder ... --genes-file my_accessions.txt     # only these accessions
    python scripts/bindingdb_local.py --folder ... --accessions P56524,Q13547,P00533

The folder should hold the unzipped files (names may carry a date, e.g. BindingDB_All_202611.tsv):
  BindingDB_All*.tsv           main table, one row per measurement          (required)
  BindingDB_rsid_eaids*.tsv    reaction-set id -> Entry ID + Assay ID        (optional, needed for assay text)
  BindingDB_Assays*.tsv        Entry ID + Assay ID -> plain-text description (optional, needed for assay text)

Proteins are matched by UniProt accession (every chain of a complex counts), never by name, so
"HDAC1" and "HDAC1/REST corepressor 3" both land on HDAC1's accession.

Writes to build_output/bindingdb/ :
  bindingdb_by_protein.csv         ligands, measurements, affinity types, sources, structures, complexes,
                                   and how many distinct assays read as inhibition / activation / unclear
  bindingdb_assays_by_protein.csv  one row per protein x assay: text, effect guess, measurement count
Column names are found by looking at the header, so small renames by BindingDB do not matter.
If something is not found the script prints the header it saw and tells you what it looked for.
Standard library only.
"""
import argparse
import collections
import csv
import glob
import html
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "build_output", "bindingdb")
AFFINITY = ["Ki", "IC50", "Kd", "EC50"]
ACC = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")
INVERSE = re.compile(r"inverse[- ]agonis\w*", re.I)
INHIB = re.compile(r"\binhibit\w*|\bantagonis\w*|\bblock\w*|\bnegative allosteric|\bdeactivat\w*", re.I)
ACTIV = re.compile(r"\bactivat\w*|\bagonis\w*|\bpotentiat\w*|\bpositive allosteric|\benhanc\w*|\bstimulat\w*|\binduc\w*", re.I)
BINDING = re.compile(r"\bbinding\b|\bdisplacement\b|\bcompetit\w*|\bradioligand\b|\baffinity\b|\bKd\b", re.I)
EFFECTS = ["inhibition", "inhibition_ic50_only", "activation", "activation_ec50_only", "binding_only", "mixed", "unclear"]
RANGE = re.compile(r"\[(\d+)-(\d+)(?:,[^\]]*)?\]")
MUTANT = re.compile(r"\b([A-Z]\d{1,4}[A-Z])\b(?=[,\]])")
REVIEWED = os.path.join(ROOT, "build_cache", "uniprot_human_reviewed.csv")


def find_file(folder, pattern, exclude=()):
    hits = [p for p in glob.glob(os.path.join(folder, "**", pattern), recursive=True)
            if not any(x.lower() in os.path.basename(p).lower() for x in exclude)
            and not p.lower().endswith((".zip", ".sdf"))]
    return sorted(hits, key=os.path.getsize, reverse=True)[0] if hits else None


def read_tsv(path):
    handle = open(path, newline="", encoding="utf-8", errors="replace")
    return handle, csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)


def pick(header, *needles, none_of=()):
    """Indices of header columns whose lower-case name contains every needle."""
    return [i for i, h in enumerate(header)
            if all(n in h.lower() for n in needles) and not any(x in h.lower() for x in none_of)]


NOISE = re.compile(r"light activation|activated sodium|sodium ortho-?vanadate|activated charcoal|activating (?:mutation|mutant)\w*"
                   r"|chelation-enhanced|enhanced fluorescence|activation loop|activation segment|constitutively activ\w*|activated (?:form|state)", re.I)
HEAD = 250        # the assay NAME and first sentence say what was measured; the protocol below is full of stray words


def effect_of(text, types):
    """Direction from the start of the assay description (its name), then from the data types. Protocol text is ignored:
    'light activation', 'activated sodium orthovanadate' and 'activating mutation' are not assay directions."""
    head = NOISE.sub(" ", text[:HEAD])
    i = bool(INHIB.search(head)) or bool(INVERSE.search(head))     # inverse agonists lower activity
    a = bool(ACTIV.search(INVERSE.sub(" ", head)))
    if i and a:
        return "mixed"
    if i:
        return "inhibition"
    if a:
        return "activation"
    if BINDING.search(head) or (types and types <= {"Ki", "Kd"}):
        return "binding_only"
    if "IC50" in types and "EC50" not in types:
        return "inhibition_ic50_only"
    if types == {"EC50"}:
        return "activation_ec50_only"
    return "unclear"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", required=True)
    ap.add_argument("--accessions", help="comma-separated UniProt accessions (default: every protein in the file)")
    ap.add_argument("--genes-file", help="text file with one accession per line")
    ap.add_argument("--min-measurements", type=int, default=1)
    ap.add_argument("--max-assays-per-protein", type=int, default=100, help="assay rows written per protein (biggest first)")
    args = ap.parse_args(argv)

    main_path = find_file(args.folder, "BindingDB_All*.tsv") or find_file(args.folder, "BindingDB*.tsv",
                                                                         exclude=("assays", "rsid", "uniprot"))
    rsid_path = find_file(args.folder, "*rsid*eaid*") or find_file(args.folder, "*rsid*")
    assay_path = find_file(args.folder, "BindingDB_Assays*")
    if not main_path:
        sys.exit("Could not find the main BindingDB_All*.tsv file inside %s" % args.folder)
    wanted = set()
    if args.accessions:
        wanted |= {a.strip().upper() for a in args.accessions.split(",") if a.strip()}
    if args.genes_file:
        with open(args.genes_file, encoding="utf-8") as h:
            wanted |= {l.strip().upper() for l in h if l.strip() and not l.startswith("#")}
    print("Main table : %s" % main_path)
    print("Assay files: %s | %s" % (rsid_path, assay_path))

    rs_to_assay, assay_text = {}, {}
    if rsid_path:
        h, r = read_tsv(rsid_path)
        next(r, None)
        for row in r:
            if len(row) >= 2:
                rs_to_assay[row[0].strip()] = row[1].strip()
        h.close()
    if assay_path:
        h, r = read_tsv(assay_path)
        next(r, None)
        for row in r:
            if len(row) >= 3:
                assay_text["%s_%s" % (row[0].strip(), row[1].strip())] = " ".join(x.strip() for x in row[2:] if x.strip())
        h.close()

    handle, reader = read_tsv(main_path)
    header = next(reader)
    c_acc = pick(header, "uniprot", "primary")
    c_rs = pick(header, "reactant_set") or pick(header, "reaction", "set")
    c_name = pick(header, "target name")
    c_src = pick(header, "curation") or pick(header, "datasource")
    c_pdb = pick(header, "pdb", "complex") or pick(header, "pdb")
    c_chains = pick(header, "number of protein chains")
    c_lig = pick(header, "monomerid") or pick(header, "ligand", "smiles")
    c_aff = {t: pick(header, t.lower() + " (nm)") for t in AFFINITY}
    problems = []
    if not c_acc:
        problems.append("no column containing 'UniProt' and 'Primary' (looked for the UniProt accession of each target chain)")
    if not c_lig:
        problems.append("no ligand id / SMILES column")
    if problems:
        print("\nHEADER SEEN (%d columns):" % len(header))
        for i, h in enumerate(header):
            print("  %3d  %s" % (i, h))
        sys.exit("\nCannot continue: " + "; ".join(problems) + ". Send me this header list and I will adjust.")
    print("Accession columns: %d | affinity columns: %s" % (len(c_acc), {t: len(v) for t, v in c_aff.items()}))

    reviewed = {r["accession"] for r in csv.DictReader(open(REVIEWED, newline="", encoding="utf-8-sig"))} \
        if os.path.exists(REVIEWED) else None
    prot = collections.defaultdict(lambda: {"lig": set(), "n": 0, "ranges": collections.Counter(), "mutants": collections.Counter(), "types": collections.Counter(),
                                            "src": collections.Counter(), "pdb": 0, "names": collections.Counter(),
                                            "partners": collections.Counter(), "single": 0, "rs": collections.defaultdict(lambda: [0, set()])})
    t0, rows = time.time(), 0
    width = len(header)
    for row in reader:
        rows += 1
        if rows % 500000 == 0:
            print("  %d rows read (%.0fs)" % (rows, time.time() - t0))
        if len(row) < width:
            row += [""] * (width - len(row))
        accs = set()
        for i in c_acc:
            accs.update(m.upper() for m in ACC.findall(row[i]))
        if reviewed:
            keep = accs & reviewed
            if keep:                      # unreviewed twins (e.g. A0A090N8Q0 = KCNH2) fold into the reviewed protein
                accs = keep
        if not accs or (wanted and not accs & wanted):
            continue
        types = {t for t, cols in c_aff.items() if any(row[i].strip() for i in cols)}
        lig = row[c_lig[0]].strip()
        rs = row[c_rs[0]].strip() if c_rs else ""
        for a in (accs & wanted if wanted else accs):
            p = prot[a]
            p["n"] += 1
            p["lig"].add(lig)
            p["types"].update(types)
            if c_src:
                p["src"][row[c_src[0]].strip()] += 1
            if c_pdb and row[c_pdb[0]].strip():
                p["pdb"] += 1
            if c_name:
                nm = row[c_name[0]].strip()
                p["names"][nm[:80]] += 1
                p["ranges"].update("%s-%s" % r for r in RANGE.findall(nm))
                p["mutants"].update(MUTANT.findall(nm))
            for other in accs - {a}:
                if reviewed is None or other in reviewed:
                    p["partners"][other] += 1
            if len(accs) == 1 and (not c_chains or row[c_chains[0]].strip() in ("", "1")):
                p["single"] += 1
            if rs:
                key = rs_to_assay.get(rs) or "rs:" + rs          # one assay can hold many measurements
                p["rs"][key][0] += 1
                p["rs"][key][1].update(types)
    handle.close()
    print("Read %d rows; %d proteins matched (%.0fs)" % (rows, len(prot), time.time() - t0))

    os.makedirs(OUT, exist_ok=True)
    prot_rows, assay_rows = [], []
    for acc, p in sorted(prot.items(), key=lambda kv: -kv[1]["n"]):
        if p["n"] < args.min_measurements:
            continue
        eff = collections.Counter()
        ordered = sorted(p["rs"].items(), key=lambda kv: -kv[1][0])
        for n_row, (key, (n, types)) in enumerate(ordered):
            text = assay_text.get(key, "")
            e = effect_of(text, types)
            eff[e] += 1
            if n_row < args.max_assays_per_protein:
                assay_rows.append({"accession": acc, "assay_key": key, "n_measurements": n,
                                   "affinity_types": ";".join(sorted(types)), "effect_guess": e, "text": html.unescape(text)[:500]})
        prot_rows.append({
            "accession": acc, "reviewed": "" if reviewed is None else ("yes" if acc in reviewed else "no (unreviewed duplicate entry?)"),
            "target_names_seen": " | ".join(n for n, _ in p["names"].most_common(3)),
            "n_measurements": p["n"], "n_single_target_measurements": p["single"], "n_ligands": len(p["lig"]),
            **{"n_" + t: p["types"][t] for t in AFFINITY}, "rows_with_structure": p["pdb"],
            "n_assays": len(p["rs"]), **{"assays_" + e: eff[e] for e in EFFECTS},
            "construct_ranges_seen": "; ".join("%s x%d" % kv for kv in p["ranges"].most_common(6)),
            "mutants_seen": "; ".join("%s x%d" % kv for kv in p["mutants"].most_common(6)),
            "top_sources": "; ".join("%s %d" % kv for kv in p["src"].most_common(3)),
            "also_listed_with": "; ".join("%s %d" % kv for kv in p["partners"].most_common(5))})
    cols = list(prot_rows[0].keys()) if prot_rows else ["accession"]
    for name, data, c in (("bindingdb_by_protein.csv", prot_rows, cols),
                          ("bindingdb_assays_by_protein.csv", assay_rows,
                           ["accession", "assay_key", "n_measurements", "affinity_types", "effect_guess", "text"])):
        with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=c)
            w.writeheader()
            w.writerows(data)
    print("Wrote %s (%d proteins, %d protein x assay rows)" % (OUT, len(prot_rows), len(assay_rows)))
    if not (rsid_path and assay_path):
        print("NOTE: assay text files not found, so effect guesses use affinity types only.")


if __name__ == "__main__":
    main()
