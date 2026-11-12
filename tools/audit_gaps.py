#!/usr/bin/env python3
"""
Where is the site thin or unreliable?  Reads the playbook, the build outputs and the website, and writes
three tables to build_output/audit/ (it changes nothing else and never uses the internet).

    py tools/audit_gaps.py
    py tools/audit_gaps.py --uniprot-tsv C:\\path\\human_reviewed.tsv.gz      # adds whole-proteome coverage

The optional --uniprot-tsv file is UniProt's table of all reviewed human proteins. Download it once with:
    Invoke-WebRequest -OutFile human_reviewed.tsv.gz -Uri "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&query=(organism_id:9606)%20AND%20(reviewed:true)&fields=accession,gene_primary,protein_name,length,ec,keyword,go_f,ft_domain,ft_topo_dom,ft_region,ft_motif,ft_repeat,ft_zn_fing,ft_dna_bind"

Writes:
  build_output/audit/gaps_by_class.csv           per function class: proteins, campaigns, direction balance, gaps
  build_output/audit/entries_missing_fields.csv  per playbook entry: what is missing (evidence, DOIs, kits ...)
  build_output/audit/unclassified_readouts.csv   readouts the website's instrument filter cannot place
(build_output/audit/vendor_sourcing_candidates.csv is written by hand, from vendor research; see GAPS.md.)
Standard library only.
"""
import argparse
import collections
import csv
import gzip
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import m4_recommend as m4          # noqa: E402  (same class rules the pipeline uses)

DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "build_output")
AUDIT = os.path.join(OUT, "audit")
FEW = 5             # the website falls back to the whole class below this many classifiable campaigns (MINC)
TSV_FEATURES = {"Domain [FT]": "Domain", "Topological domain": "Topological domain", "Region": "Region",
                "Motif": "Motif", "Repeat": "Repeat", "Zinc finger": "Zinc finger", "DNA binding": "DNA binding"}


def read(path):
    try:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            return list(csv.DictReader(handle))
    except OSError:
        return []


def write(path, rows, columns):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        w = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def split(cell):
    return [x.strip() for x in re.split(r"[;|]", cell or "") if x.strip()]


def is_cell(scope):
    return "cell-based" in (scope or "") and "free" not in (scope or "")


# ------------------------------------------------------------------ the website's instrument keywords

def instrument_rules():
    """The INSTR list from docs/index.html, so this check always matches what the site does."""
    with open(os.path.join(ROOT, "docs", "index.html"), encoding="utf-8") as handle:
        page = handle.read()
    block = page[page.index("const INSTR=["):page.index("const instOf")]
    return [(k, re.compile(rx, re.I)) for k, rx in re.findall(r'\["(\w+)",/((?:\\/|[^/])+)/i\]', block)]


def instruments_for(text, rules):
    return [k for k, rx in rules if rx.search(text)]


# ------------------------------------------------------------------ whole proteome (optional)

def proteome_classes(tsv_path, class_rows):
    """accession -> list of function classes, using the pipeline's own rules on UniProt's annotation table."""
    out = {}
    with gzip.open(tsv_path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            feats = []
            for col, ftype in TSV_FEATURES.items():
                for m in re.finditer(r"(\d+)\.\.(\d+);\s*/note=\"([^\"]*)\"", row.get(col, "") or ""):
                    feats.append({"type": ftype, "start": int(m.group(1)), "end": int(m.group(2)), "description": m.group(3)})
            info = {"ecs": [e.strip() for e in (row.get("EC number") or "").split(";") if e.strip()],
                    "keywords": [k.strip() for k in (row.get("Keywords") or "").split(";") if k.strip()],
                    "go": [re.sub(r"\s*\[GO:\d+\]", "", g).strip() for g in (row.get("Gene Ontology (molecular function)") or "").split(";") if g.strip()],
                    "features": feats}
            found = m4.classify(info, class_rows)
            out[row["Entry"]] = {"gene": row.get("Gene Names (primary)", ""), "name": row.get("Protein names", ""),
                                 "classes": [c for c in found if c != "generic"]}
    return out


# ------------------------------------------------------------------ main

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--uniprot-tsv", help="UniProt table of reviewed human proteins (.tsv.gz) for whole-proteome coverage")
    args = ap.parse_args(argv)

    assays = read(os.path.join(DATA, "assays.csv"))
    kits = read(os.path.join(DATA, "kits.csv"))
    class_rows = m4.read_table(os.path.join(DATA, "function_classes.csv"))
    usage = read(os.path.join(DATA, "usage_stats.csv"))
    by_class = {r["class_id"]: r for r in read(os.path.join(OUT, "effect_by_class.csv"))}
    gaps = {r["class_id"]: r for r in read(os.path.join(OUT, "function_gaps.csv"))}
    genes = read(os.path.join(OUT, "effect_by_gene.csv"))
    suggestions = read(os.path.join(OUT, "readout_suggestions.csv"))
    rules = instrument_rules()

    # --- per entry ------------------------------------------------------------------------------------------
    kits_by_entry = collections.defaultdict(list)
    for k in kits:
        for e in split(k["applies_to_entries"]):
            kits_by_entry[e].append(k)
    usage_by_entry = collections.defaultdict(list)
    for u in usage:
        usage_by_entry[u["entry_id"]].append(u)
    groups = collections.defaultdict(list)          # entries that can only ever show the same headline share
    for a in assays:
        ev = tuple(sorted(split(a.get("evidence_labels") or a.get("evidence_function"))))
        if ev:
            for c in split(a["classes"]):
                groups[(c, is_cell(a["scope"]), ev)].append(a["entry_id"])
    twins = collections.defaultdict(set)
    for members in groups.values():
        for e in members:
            twins[e].update(m for m in members if m != e)

    entry_rows, unclassified = [], []
    for a in assays:
        e = a["entry_id"]
        ks = kits_by_entry.get(e, [])
        rows = usage_by_entry.get(e, [])
        tops = [u for u in rows if not u["readout"]]
        genes_ok = {u["scope_id"] for u in tops if u["scope_type"] == "gene" and int(u["n_campaigns_classifiable"] or 0) >= FEW}
        readouts = [r for r in split(a["readouts"]) if not r.lower().startswith("same")]
        unknown = [r for r in readouts if not instruments_for(r, rules)]
        for r in unknown:
            unclassified.append({"entry_id": e, "readout": r, "source": "playbook readout list"})
        missing = {
            "no_evidence_labels": not split(a.get("evidence_labels") or a.get("evidence_function")),
            "no_sources_doi": not split(a.get("sources_doi")),
            "draft": "draft" in (a.get("status") or "").lower(),
            "no_kit": not ks,
            "only_unverified_kits": bool(ks) and all(k.get("status") == "unverified" for k in ks),
            "no_plate_formats": not split(a.get("plate_formats")),
            "no_learn_more_url": not a.get("learn_more_url"),
            "no_usage_rows": not rows,
        }
        entry_rows.append(dict(
            {"entry_id": e, "assay_name": a["assay_name"], "classes": a["classes"], "role": a["role"],
             "status": a["status"], "n_kits": len(ks), "n_usage_rows": len(rows),
             "n_proteins_with_5plus_campaigns": len(genes_ok),
             "same_headline_share_as": "; ".join(sorted(twins.get(e, ()))),
             "readouts_instrument_filter_cannot_place": "; ".join(unknown),
             "n_missing": sum(missing.values())},
            **{k: "yes" if v else "" for k, v in missing.items()}))
    entry_rows.sort(key=lambda r: (-r["n_missing"], r["entry_id"]))
    write(os.path.join(AUDIT, "entries_missing_fields.csv"), entry_rows,
          ["entry_id", "assay_name", "classes", "role", "status", "n_missing", "no_evidence_labels", "no_sources_doi",
           "draft", "no_kit", "only_unverified_kits", "no_plate_formats", "no_learn_more_url", "no_usage_rows",
           "n_kits", "n_usage_rows", "n_proteins_with_5plus_campaigns", "same_headline_share_as",
           "readouts_instrument_filter_cannot_place"])

    # technologies seen in the data that the instrument filter cannot place
    seen = collections.Counter()
    for s in suggestions:
        seen[s["technology"]] += int(s["n_campaigns"] or 0)
    for tech, n in seen.most_common():
        if not instruments_for(tech, rules):
            unclassified.append({"entry_id": "", "readout": tech, "source": "data (readout_suggestions.csv), %d campaign-pairs" % n})
    write(os.path.join(AUDIT, "unclassified_readouts.csv"), unclassified, ["entry_id", "readout", "source"])

    # --- per class ------------------------------------------------------------------------------------------
    entries_by_class = collections.defaultdict(list)
    for a in assays:
        for c in split(a["classes"]):
            entries_by_class[c].append(a)
    genes_by_class = collections.defaultdict(list)
    for g in genes:
        for c in split(g["classes"]):
            genes_by_class[c].append(g)
    not_covered = collections.Counter(s["class_id"] for s in suggestions if s["suggestion"] != "covered")
    proteome = proteome_classes(args.uniprot_tsv, class_rows) if args.uniprot_tsv else {}
    with_data = {g["accession"] for g in genes if int(g["n_campaigns"] or 0) > 0}
    proteome_by_class = collections.defaultdict(set)
    for acc, p in proteome.items():
        for c in p["classes"] or ["(no class)"]:
            proteome_by_class[c].add(acc)

    class_out = []
    for row in class_rows + [{"class_id": "(no class)", "label": "No function class matched"}]:
        c = row["class_id"]
        if c == "generic":
            continue
        e = by_class.get(c, {})
        n_c = int(e.get("n_campaigns") or 0)
        act = int(e.get("activation") or 0) + int(e.get("mixed") or 0)
        inh = n_c - int(e.get("activation") or 0)
        gs = genes_by_class.get(c, [])
        usable = [int(g["n_campaigns"] or 0) - int(g["no_function_label"] or 0) for g in gs]
        ents = entries_by_class.get(c, [])
        members = proteome_by_class.get(c, set())
        class_out.append({
            "class_id": c, "label": row.get("label", ""),
            "n_playbook_entries": len(ents),
            "n_entries_reviewed": sum(1 for a in ents if "draft" not in a["status"].lower()),
            "n_entries_with_kit": sum(1 for a in ents if kits_by_entry.get(a["entry_id"])),
            "n_entries_without_evidence_labels": sum(1 for a in ents if not split(a.get("evidence_labels") or a.get("evidence_function"))),
            "n_proteins_in_proteome": len(members) if proteome else "",
            "n_proteins_in_proteome_with_campaigns": len(members & with_data) if proteome else "",
            "n_proteins_in_proteome_zero_campaigns": len(members - with_data) if proteome else "",
            "n_proteins_summarized": len(gs),
            "n_proteins_lt5_usable_campaigns": sum(1 for u in usable if u < FEW),
            "n_campaigns": n_c,
            "campaigns_for_inhibition_view": inh, "campaigns_for_activation_view": act,
            "activation_share": round(act / n_c, 3) if n_c else "",
            "pct_no_function_label": round(int(gaps.get(c, {}).get("no_function_label") or 0) / n_c, 3) if n_c else "",
            "pct_no_technology_label": round(int(gaps.get(c, {}).get("no_technology_label") or 0) / n_c, 3) if n_c else "",
            "readout_suggestions_not_in_playbook": not_covered.get(c, 0),
        })
    write(os.path.join(AUDIT, "gaps_by_class.csv"), class_out, list(class_out[0].keys()))

    # --- printed summary (used for GAPS.md) -------------------------------------------------------------------
    print("Playbook: %d entries; missing evidence_labels %d, no DOI %d, draft %d, no kit %d, no plate formats %d, no usage rows %d"
          % (len(assays), *[sum(1 for r in entry_rows if r[k]) for k in
                             ("no_evidence_labels", "no_sources_doi", "draft", "no_kit", "no_plate_formats", "no_usage_rows")]))
    print("Entries that always show the same headline share as a sibling: %d" % sum(1 for r in entry_rows if r["same_headline_share_as"]))
    print("Playbook readouts the instrument filter cannot place: %d (in %d entries)"
          % (sum(1 for u in unclassified if u["entry_id"]), len({u["entry_id"] for u in unclassified if u["entry_id"]})))
    if proteome:
        zero = [a for a in proteome if a not in with_data]
        print("Proteome: %d reviewed human proteins; %d with >=1 campaign in our data; %d with none"
              % (len(proteome), len(set(proteome) & with_data), len(zero)))
        cz = collections.Counter(c for a in zero for c in (proteome[a]["classes"] or ["(no class)"]))
        print("  zero-campaign proteins by class: " + ", ".join("%s %d" % kv for kv in cz.most_common()))
        write(os.path.join(AUDIT, "proteome_zero_campaigns.csv"),
              [{"accession": a, "gene": proteome[a]["gene"], "classes": "; ".join(proteome[a]["classes"]),
                "name": proteome[a]["name"][:120]} for a in sorted(zero)], ["accession", "gene", "classes", "name"])
    print("Wrote %s" % os.path.relpath(AUDIT, ROOT))


if __name__ == "__main__":
    main()
