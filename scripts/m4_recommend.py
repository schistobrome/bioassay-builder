"""
M4: assay recommendations from UniProt + a curated assay playbook.

You give a protein (UniProt accession or gene symbol) and, optionally, the
residues your construct covers. The script
  1. reads the protein's UniProt entry (domains with residue ranges, EC numbers,
     keywords, GO terms, PDB structures),
  2. decides which FUNCTION CLASSES the protein belongs to (kinase, helicase,
     nuclease, protease, phosphatase, receptor ectodomain, plus a generic set),
  3. checks how much of each class's domain your construct covers, and drops the
     assays whose domain is outside the construct,
  4. prints the playbook assays that remain, grouped as primary / confirmation /
     counterscreen, with readouts, typical conditions and a reagent checklist,
  5. if m2_data exists, adds how often the public data (ChEMBL/PubChem) used each
     assay type, as supporting evidence only.

The curated tables are plain CSV files that YOU can edit in Excel:
    function_classes.csv   which proteins count as which class
    assay_playbook.csv     the assays for each class (first created by this script
                           the reviewed version ships with this script)
They are created on the first run. Use --reset-playbook to replace them with the
newest starter versions (the old ones are kept).

How to run (with the Python in .venv, in the AssayBuilder folder):
    .\\.venv\\Scripts\\python.exe m4_recommend.py WRN
    .\\.venv\\Scripts\\python.exe m4_recommend.py WRN --start 517 --end 1093
    .\\.venv\\Scripts\\python.exe m4_recommend.py EGFR --start 696 --end 1022
    .\\.venv\\Scripts\\python.exe m4_recommend.py P00533 --start 25 --end 645
Needs internet only the first time for a protein (answers are saved in m4_data).
It only READS public data. Uses only Python's built-in libraries.
"""

import argparse
import collections
import csv
import json
import os
import re
import sys
import urllib.parse
import urllib.request
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
CLASSES_CSV = os.path.join(HERE, "function_classes.csv")
KITS_CSV = os.path.join(HERE, "commercial_kits.csv")   # vendor kits linked to playbook entries
PLAYBOOK_CSV = os.path.join(HERE, "assay_playbook.csv")   # curated; ship it next to this script
CACHE = os.path.join(HERE, "m4_data")
M1 = os.path.join(HERE, "m1_data")
M2 = os.path.join(HERE, "m2_data")
HEADERS = {"User-Agent": "assay-planner-m4/0.1 (personal side project)", "Accept": "application/json"}
SHORTCUTS = {"WRN": "Q14191", "EGFR": "P00533"}
FULL, PARTIAL = 0.9, 0.0      # coverage of a class's domain: >= FULL is covered, > PARTIAL is partial
ROLE_ORDER = ["primary", "confirmation", "counterscreen"]
CATEGORY_LABEL = {
    "biochemical_hts": "Established biochemical HTS assay",
    "cell_functional": "Cell-based functional screening assay",
    "biophysical_binding": "Biophysical binding / target-engagement assay",
    "orthogonal_confirmation": "Orthogonal confirmation assay",
    "mechanistic_counterscreen": "Mechanistic / selectivity counterscreen",
}
LEVEL_LABEL = {"very_low": "very low", "low": "low", "moderate": "moderate", "high": "high", "very_high": "very high"}
ANCHOR_TYPES = {"Domain", "Region", "Repeat", "Motif", "Zinc finger", "Topological domain", "DNA binding"}

LOG = []


def say(text=""):
    print(text)
    LOG.append(text)


# ---------------------------------------------------------------- starter tables

CLASS_COLUMNS = ["class_id", "label", "ec_prefixes", "keywords", "go_terms", "domain_patterns",
                 "needs_domain", "modality", "excludes", "notes"]
NUC_EC = "3.1.11|3.1.13|3.1.14|3.1.15|3.1.16|3.1.21|3.1.22|3.1.25|3.1.26|3.1.27|3.1.30|3.1.31|3.1.-"
CLASSES = [
    # class_id, label, EC prefixes, UniProt keywords, GO function terms, domain-name patterns,
    # needs_domain, modality flag, classes that cancel this one, notes   (lists separated by | or ;)
    ("protein_kinase", "Protein kinase", "2.7.10|2.7.11|2.7.12|2.7.13", "kinase",
     "protein kinase activity", "protein kinase", "no", "", "", ""),
    ("dna_helicase", "ATP-dependent helicase", "3.6.4.12|3.6.4.13", "helicase",
     "helicase activity", "helicase", "no", "", "", ""),
    ("nuclease", "Nuclease (exo- / endonuclease)", NUC_EC, "exonuclease|endonuclease|nuclease",
     "exonuclease activity|endonuclease activity|nuclease activity", "exonuclease|endonuclease|nuclease",
     "no", "", "", ""),
    ("protease", "Protease / peptidase", "3.4.", "protease", "peptidase activity",
     "peptidase|protease|caspase", "no", "", "", ""),
    ("protein_phosphatase", "Protein phosphatase", "3.1.3.16|3.1.3.48", "protein phosphatase",
     "protein tyrosine phosphatase activity|protein serine/threonine phosphatase activity",
     "phosphatase", "no", "", "", ""),
    ("receptor_ectodomain", "Receptor ligand-binding ectodomain", "", "receptor", "receptor activity",
     "extracellular", "yes", "", "gpcr|ion_channel|abc_transporter|solute_carrier",
     "Only counted when UniProt annotates an extracellular region (Topological domain)."),
    ("gpcr", "G-protein-coupled receptor", "", "G-protein coupled receptor",
     "G protein-coupled receptor activity|G-protein coupled receptor activity", "", "no", "", "", ""),
    ("ion_channel", "Ion channel", "",
     "Ion channel|Voltage-gated channel|Ligand-gated ion channel|Potassium channel|Sodium channel|Calcium channel|Chloride channel",
     "ion channel activity|channel activity", "", "no", "", "", ""),
    ("abc_transporter", "ABC transporter", "", "ABC transporter|Multidrug resistance",
     "ABC type transporter activity|ATPase coupled transmembrane transporter activity",
     "abc transporter|abc transmembrane", "no", "", "", ""),
    ("solute_carrier", "Solute carrier / secondary transporter", "",
     "Symport|Antiport|Sugar transport|Amino acid transport|Neurotransmitter transport",
     "symporter activity|antiporter activity|secondary active transmembrane transporter activity", "",
     "no", "", "", ""),
    ("epigenetic_writer", "Histone methyltransferase (writer)",
     "2.1.1.43|2.1.1.319|2.1.1.320|2.1.1.321|2.1.1.322|2.1.1.367", "",
     "histone methyltransferase activity|histone arginine methyltransferase activity",
     "set|dot1|prmt|protein arginine n methyltransferase", "no", "", "", ""),
    ("epigenetic_demethylase", "Histone demethylase (eraser)",
     "1.14.11.27|1.14.11.65|1.14.11.66|1.14.11.67|1.14.11.68|1.14.11.69|1.14.99.66", "",
     "histone demethylase activity", "jmjc", "no", "", "", ""),
    ("epigenetic_deacetylase", "Histone deacetylase / sirtuin (eraser)", "3.5.1.98|2.3.1.286", "",
     "histone deacetylase activity|NAD dependent protein deacetylase activity",
     "histone deacetylase|sirtuin", "no", "", "", ""),
    ("epigenetic_reader", "Bromodomain reader", "", "Bromodomain", "", "bromo|bromodomain", "no", "", "", ""),
    ("ppi_domain", "Protein-protein interaction domain", "", "", "",
     "sh2|sh3|pdz|ww|wd|ank|tpr|btb|poz|bh3", "no", "", "",
     "Broad on purpose: any protein with one of these interaction domains gets the PPI assays."),
    ("e3_ligase", "E3 ubiquitin ligase", "2.3.2.27|2.3.2.26|2.3.2.31", "",
     "ubiquitin protein ligase activity|ubiquitin-protein transferase activity", "hect|rbr|u box", "no", "", "", ""),
    ("dub", "Deubiquitinase", "3.4.19.12", "",
     "ubiquitin specific protease activity|cysteine type deubiquitinase activity|deubiquitinase activity",
     "usp|uch|josephin|otu|mindy", "no", "", "", ""),
    ("degrader_modality", "Targeted protein degradation project", "", "", "", "", "no", "degrader", "",
     "Only included when you pass --modality degrader."),
    ("nuclear_receptor", "Nuclear receptor", "", "",
     "nuclear receptor activity|steroid hormone receptor activity", "nuclear receptor|nr c4 type|nr lbd", "no", "", "", ""),
    ("topoisomerase", "DNA topoisomerase", "5.6.2.1|5.6.2.2|5.99.1.2|5.99.1.3", "Topoisomerase",
     "DNA topoisomerase activity", "topoisomerase|toprim", "no", "", "", ""),
    ("rna_binding", "RNA-binding protein", "", "RNA-binding", "RNA binding",
     "rrm|kh|dsrm|pumilio", "no", "rna_target", "",
     "UniProt does not annotate the RNA itself; for RNA-structure targets pass --modality rna_target."),
    ("autophagy", "Autophagy / lysosomal flux", "", "Autophagy", "", "", "no", "autophagy", "",
     "Also available for any protein with --modality autophagy."),
    ("generic", "Any protein (biophysical and interference checks)", "", "", "", "", "no", "", "",
     "Always applies."),
]


def write_table(path, columns, rows):
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows(rows)


def ensure_tables(reset):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if reset and os.path.exists(CLASSES_CSV):
        os.replace(CLASSES_CSV, CLASSES_CSV.replace(".csv", "_old_%s.csv" % stamp))
        say("Old function_classes.csv kept as a dated copy.")
    if os.path.exists(CLASSES_CSV):
        with open(CLASSES_CSV, encoding="utf-8-sig") as handle:
            header = handle.readline()
        if "modality" not in header:
            os.replace(CLASSES_CSV, CLASSES_CSV.replace(".csv", "_old_%s.csv" % stamp))
            say("function_classes.csv was the older layout; replaced with the new one (old copy kept).")
    if not os.path.exists(CLASSES_CSV):
        write_table(CLASSES_CSV, CLASS_COLUMNS, CLASSES)
        say("Created the starter table: %s" % CLASSES_CSV)
    if not os.path.exists(PLAYBOOK_CSV):
        raise SystemExit("assay_playbook.csv was not found next to this script (%s). Put the playbook file there." % PLAYBOOK_CSV)


def read_table(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


# ---------------------------------------------------------------- UniProt

def http_get_text(url, timeout=60):
    request = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def resolve_accession(query):
    q = query.strip()
    if q.upper() in SHORTCUTS:
        return SHORTCUTS[q.upper()]
    if re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", q.upper()):
        return q.upper()
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": "gene_exact:%s AND organism_id:9606 AND reviewed:true" % q, "format": "json", "size": 1, "fields": "accession"})
    data = json.loads(http_get_text(url))
    if not data.get("results"):
        raise SystemExit("No human reviewed UniProt entry found for '%s'. Try the accession (like Q14191)." % q)
    return data["results"][0]["primaryAccession"]


def load_entry(accession):
    for folder in (M1, CACHE):
        path = os.path.join(folder, "uniprot_%s.json" % accession)
        if os.path.exists(path):
            with open(path, encoding="utf-8") as handle:
                return json.load(handle), path
    os.makedirs(CACHE, exist_ok=True)
    text = http_get_text("https://rest.uniprot.org/uniprotkb/%s.json" % accession)
    path = os.path.join(CACHE, "uniprot_%s.json" % accession)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return json.loads(text), path


def find_values(obj, key, out=None):
    out = [] if out is None else out
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k == key and isinstance(v, list):
                out.extend(x.get("value", "") for x in v if isinstance(x, dict))
            else:
                find_values(v, key, out)
    elif isinstance(obj, list):
        for item in obj:
            find_values(item, key, out)
    return out


def norm(text):
    return re.sub(r"[^a-z0-9]+", " ", str(text).lower()).strip()


def summarize_entry(entry):
    gene = ""
    if entry.get("genes"):
        gene = (entry["genes"][0].get("geneName") or {}).get("value", "")
    length = int((entry.get("sequence") or {}).get("length") or 0)
    ecs = sorted(set(v for v in find_values(entry.get("proteinDescription", {}), "ecNumbers") if v))
    keywords = [k.get("name", "") for k in entry.get("keywords", [])]
    go = []
    pdb = []
    for ref in entry.get("uniProtKBCrossReferences", []):
        props = {p.get("key"): p.get("value") for p in ref.get("properties", [])}
        if ref.get("database") == "GO" and props.get("GoTerm", "").startswith("F:"):
            go.append(props["GoTerm"][2:])
        if ref.get("database") == "PDB":
            for m in re.finditer(r"=(\d+)-(\d+)", props.get("Chains", "")):
                pdb.append((ref.get("id"), int(m.group(1)), int(m.group(2))))
    features = []
    for f in entry.get("features", []):
        if f.get("type") in ANCHOR_TYPES:
            loc = f.get("location", {})
            try:
                features.append({"type": f["type"], "description": f.get("description", ""),
                                 "start": int(loc["start"]["value"]), "end": int(loc["end"]["value"])})
            except (KeyError, TypeError, ValueError):
                pass
    return {"gene": gene, "length": length, "ecs": ecs, "keywords": keywords, "go": go,
            "pdb": pdb, "features": features,
            "name": ((entry.get("proteinDescription") or {}).get("recommendedName") or {}).get(
                "fullName", {}).get("value", "")}


# ---------------------------------------------------------------- classes and coverage

def split(cell):
    return [x.strip() for x in re.split(r"[|;]", cell or "") if x.strip()]


def word_in(pattern, text):
    return bool(re.search(r"(?<![a-z0-9])%s(?![a-z0-9])" % re.escape(pattern), text))


def classify(info, classes, modalities=()):
    """Return {class_id: {label, why, anchors}} for classes this protein belongs to."""
    found = collections.OrderedDict()
    for row in classes:
        cid = row["class_id"]
        if cid == "generic":
            found[cid] = {"label": row["label"], "why": ["applies to every protein"], "anchors": []}
            continue
        why = []
        if row.get("modality") and row["modality"] in modalities:
            why.append("requested with --modality %s" % row["modality"])
        for prefix in split(row["ec_prefixes"]):
            hit = [e for e in info["ecs"] if e.startswith(prefix)]
            if hit:
                why.append("EC %s" % hit[0])
                break
        for kw in split(row["keywords"]):
            if any(norm(kw) == norm(k) for k in info["keywords"]):
                why.append("UniProt keyword '%s'" % kw)
                break
        for term in split(row["go_terms"]):
            hit = [g for g in info["go"] if norm(term) in norm(g)]
            if hit:
                why.append("GO '%s'" % hit[0])
                break
        anchors = []
        patterns = [norm(p) for p in split(row["domain_patterns"])]
        for feat in info["features"]:
            if any(p and word_in(p, norm(feat["description"])) for p in patterns):
                if row["needs_domain"] == "yes" and feat["type"] != "Topological domain":
                    continue
                anchors.append(feat)
        if anchors:
            why.append("domain '%s'" % anchors[0]["description"])
        if row["needs_domain"] == "yes" and not anchors:
            continue
        if why:
            found[cid] = {"label": row["label"], "why": why, "anchors": anchors}
    for row in classes:                      # some classes cancel others (e.g. a GPCR is not a plain ectodomain receptor)
        if row["class_id"] in found:
            if any(ex in found for ex in split(row.get("excludes"))):
                del found[row["class_id"]]
    return found


def coverage(anchors, start, end):
    """Fraction of the class's domain residues inside the construct (union of its domains)."""
    if not anchors:
        return None
    spans = sorted((a["start"], a["end"]) for a in anchors)
    merged = []
    for s, e in spans:
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    total = sum(e - s + 1 for s, e in merged)
    inside = sum(max(0, min(e, end) - max(s, start) + 1) for s, e in merged)
    return inside / float(total), merged


def kits_for(kits, entry_id, gene):
    """Kits whose entry list includes this assay; target-specific kits only when the gene symbol matches."""
    out = []
    for k in kits:
        if entry_id not in split(k.get("applies_to_entries", "")):
            continue
        scope = k.get("target_scope") or ""
        if scope.lower().startswith(("pan-class", "many")):
            out.append(k)
        elif gene and re.search(r"\b%s\b" % re.escape(gene), scope, re.IGNORECASE):
            out.append(k)
    return out


def clean_scope(row):
    return "cell-based" if "cell" in (row.get("scope") or "").lower() and "free" not in (row.get("scope") or "").lower() else "cell-free"


def main_role(row):
    text = (row.get("role") or "").lower()
    for role in ROLE_ORDER:
        if role in text:
            return role
    return "confirmation"


def priority_number(row):
    """'P2/P3' -> 2, 'P1 (binding campaigns)' -> 1; unknown -> 9."""
    match = re.search(r"P\s*(\d)", row.get("development_priority") or "")
    return int(match.group(1)) if match else 9


# ---------------------------------------------------------------- public-data evidence

def evidence_counts(gene, entries):
    path = os.path.join(M2, "%s_labeled.csv" % gene)
    if not gene or not os.path.exists(path):
        return {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        rows = [r for r in csv.DictReader(handle) if r.get("in_v1_scope") == "yes"]
    result = {}
    for entry in entries:
        functions = split(entry.get("evidence_labels") or entry.get("evidence_function"))
        if not functions:
            continue
        camps, readouts = set(), collections.defaultdict(set)
        for r in rows:
            labels = [x.strip() for x in r.get("function_labels", "").split("|")]
            if any(f in labels for f in functions):
                camps.add(r["campaign"])
                for tech in [x.strip() for x in r.get("technology_labels", "").split("|") if x.strip()]:
                    readouts[tech].add(r["campaign"])
        result[entry["entry_id"]] = (len(camps), sorted(((t, len(c)) for t, c in readouts.items()),
                                                        key=lambda kv: -kv[1]))
    return result


# ---------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("protein", help="UniProt accession (e.g. Q14191) or human gene symbol (e.g. WRN)")
    parser.add_argument("--start", type=int, default=None, help="first residue of your construct")
    parser.add_argument("--end", type=int, default=None, help="last residue of your construct")
    parser.add_argument("--effect", default="inhibition", choices=["inhibition", "activation"],
                        help="desired effect (default inhibition)")
    parser.add_argument("--reset-classes", dest="reset_playbook", action="store_true",
                        help="replace function_classes.csv with the newest starter version (old one is kept)")
    parser.add_argument("--modality", nargs="*", default=[],
                        help="project types to switch on: degrader, autophagy, rna_target")
    parser.add_argument("--category", nargs="*", default=[],
                        help="show only these categories: biochemical_hts cell_functional biophysical_binding "
                             "orthogonal_confirmation mechanistic_counterscreen")
    parser.add_argument("--cell-free-only", action="store_true",
                        help="hide cell-based assays (by default they are shown and labeled)")
    parser.add_argument("--include-cell-based", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--no-evidence", action="store_true", help="skip the public-data counts")
    args = parser.parse_args()

    ensure_tables(args.reset_playbook)
    classes = read_table(CLASSES_CSV)
    playbook = read_table(PLAYBOOK_CSV)
    kits = read_table(KITS_CSV) if os.path.exists(KITS_CSV) else []
    need = {"entry_id", "assay_name", "classes", "role", "readouts", "measures", "conditions", "reagents"}
    if not playbook or not need <= set(playbook[0].keys()):
        raise SystemExit("assay_playbook.csv is missing columns: %s" % ", ".join(sorted(need - set(playbook[0].keys() if playbook else []))))

    accession = resolve_accession(args.protein)
    entry, source = load_entry(accession)
    info = summarize_entry(entry)
    start = args.start or 1
    end = args.end or info["length"]
    if start > end or (info["length"] and (start < 1 or end > info["length"])):
        raise SystemExit("Construct %d-%d does not fit the protein (1-%d)." % (start, end, info["length"]))

    say("M4 assay recommendation (entries marked 'draft by Claude' have not been expert-reviewed)")
    say("Run at: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    say("")
    say("Protein: %s (%s), %s, %d residues" % (info["gene"] or "?", accession, info["name"] or "?", info["length"]))
    say("Construct: residues %d-%d%s   Desired effect: %s" % (
        start, end, " (full length)" if start == 1 and end == info["length"] else "", args.effect))
    say("UniProt gave: EC %s | %d keywords | %d GO functions | %d PDB residue ranges"
        % (", ".join(info["ecs"]) or "none", len(info["keywords"]), len(info["go"]), len(info["pdb"])))

    found = classify(info, classes, args.modality)
    say("")
    say("Function classes found for this protein:")
    status = {}
    for cid, c in found.items():
        cov = coverage(c["anchors"], start, end)
        if cid == "generic":
            status[cid] = ("applies", None)
            say("  - %-28s %s" % (c["label"], "; ".join(c["why"])))
            continue
        if cov is None:
            status[cid] = ("applies", None)
            say("  - %-28s %s; no domain range in UniProt, so construct coverage cannot be checked" % (
                c["label"], "; ".join(c["why"])))
            continue
        fraction, spans = cov
        ranges = ", ".join("%d-%d" % tuple(s) for s in spans)
        if fraction >= FULL:
            status[cid] = ("applies", fraction)
            note = "domain %s fully inside your construct" % ranges
        elif fraction > PARTIAL:
            status[cid] = ("partial", fraction)
            note = "WARNING: only %d%% of domain %s is inside your construct" % (round(100 * fraction), ranges)
        else:
            status[cid] = ("outside", 0.0)
            note = "domain %s is OUTSIDE your construct, so these assays are dropped" % ranges
        say("  - %-28s %s; %s" % (c["label"], "; ".join(c["why"]), note))

    gene_evidence = {} if args.no_evidence else evidence_counts(info["gene"], playbook)
    chosen, dropped, out_of_scope = [], [], []
    for row in playbook:
        entry_classes = [c for c in split(row["classes"]) if c in status]
        if not entry_classes:
            continue
        states = [status[c][0] for c in entry_classes]
        if "applies" in states:
            level = "applies"
        elif "partial" in states:
            level = "partial"
        else:
            dropped.append(row)
            continue
        if args.category and row.get("assay_category") not in args.category:
            continue
        if clean_scope(row) == "cell-based" and args.cell_free_only:
            out_of_scope.append(row)
            continue
        chosen.append((row, level, entry_classes))

    counts = collections.Counter(x[0].get("assay_category") for x in chosen)
    say("")
    say("Assays shown by category: " + "; ".join(
        "%s %d" % (CATEGORY_LABEL[k], counts.get(k, 0)) for k in CATEGORY_LABEL))
    chosen.sort(key=lambda x: (ROLE_ORDER.index(main_role(x[0])), priority_number(x[0])))
    for role in ROLE_ORDER:
        group = [x for x in chosen if main_role(x[0]) == role]
        say("")
        say("=" * 78)
        say("%s ASSAYS (%d)" % (role.upper(), len(group)))
        say("=" * 78)
        if not group:
            say("")
            say("(none for this protein with the current filters%s)" % (
                "; cell-based assays are hidden by --cell-free-only" if args.cell_free_only else ""))
        for n, (row, level, entry_classes) in enumerate(group, 1):
            say("")
            say("%d. %s   [%s]%s" % (n, row["assay_name"], row["entry_id"],
                                      "   ** domain only partly in construct **" if level == "partial" else ""))
            say("   Category : %s%s" % (CATEGORY_LABEL.get(row.get("assay_category"), row.get("assay_category") or "-"),
                                          " (cell-based)" if clean_scope(row) == "cell-based" else ""))
            say("   Develop priority: %s   |   HTS suitability: %s   |   Biological relevance: %s" % (
                row.get("development_priority") or "-",
                LEVEL_LABEL.get(row.get("screening_suitability"), row.get("screening_suitability") or "-"),
                LEVEL_LABEL.get(row.get("biological_relevance"), row.get("biological_relevance") or "-")))
            if row.get("rating_notes"):
                say("   Rating notes: %s" % row["rating_notes"])
            say("   Measures : %s" % row["measures"])
            say("   Rationale: %s  (for: %s)" % (row.get("rationale") or row.get("why", ""),
                                                 ", ".join(found[c]["label"] for c in entry_classes)))
            say("   Readouts : %s" % "; ".join(split(row["readouts"])) if "validated" not in row["readouts"].lower()
                else "   Readouts : %s" % row["readouts"])
            if row["conditions"]:
                say("   Typical conditions: %s" % row["conditions"])
            say("   Reagents : %s" % row["reagents"])
            ev = gene_evidence.get(row["entry_id"])
            if ev:
                top = ", ".join("%s %d" % kv for kv in ev[1][:4]) or "no readout stated in the records"
                say("   Public data (supporting only): %d independent campaigns labeled with this function; "
                    "readouts seen: %s" % (ev[0], top))
            my_kits = kits_for(kits, row["entry_id"], info["gene"])
            if my_kits:
                say("   Commercial options (listings seen on the web; confirm on the vendor datasheet):")
                for k in my_kits:
                    say("     - %s: %s%s [%s]%s" % (
                        k["vendor"], k["product"],
                        " (%s)" % k["product_code"] if k.get("product_code") and k["product_code"] != "see vendor" else "",
                        k.get("target_scope") or "-",
                        "  ** UNVERIFIED: user-reported, not yet checked on the web **" if k.get("status") == "unverified" else ""))
                    if k.get("url"):
                        say("       %s" % k["url"])
            if row["entry_id"] == "G4" and info["pdb"]:
                anchors = [a for cid in status if status[cid][0] != "outside" and cid in found
                           for a in (found[cid].get("anchors") or [])]
                spans = [(a["start"], a["end"]) for a in anchors if min(a["end"], end) >= max(a["start"], start)]
                allp = {p[0] for p in info["pdb"]}
                if spans:
                    def hits(pid):
                        for pid2, s_, e_ in info["pdb"]:
                            if pid2 != pid:
                                continue
                            for a_s, a_e in spans:
                                lo, hi = max(a_s, start), min(a_e, end)
                                if max(0, min(e_, hi) - max(s_, lo) + 1) >= 0.5 * (hi - lo + 1):
                                    return True
                        return False
                    covering = sorted(p for p in allp if hits(p))
                    what = "cover at least half of an annotated domain in your construct"
                else:
                    span = end - start + 1
                    covering = sorted({pid for pid, s_, e_ in info["pdb"]
                                       if max(0, min(e_, end) - max(s_, start) + 1) >= 0.5 * span})
                    what = "cover at least half of your construct"
                say("   Structures: UniProt lists %d PDB entries; %d %s%s"
                    % (len(allp), len(covering), what,
                       " (e.g. %s)" % ", ".join(covering[:4]) if covering else ""))
            dois = split(row.get("sources_doi", ""))
            if dois:
                say("   Sources  : %s" % "; ".join(d if d.startswith("http") else "https://doi.org/" + d for d in dois))
            if row.get("plate_formats"):
                say("   Plate formats: %s" % "; ".join(split(row["plate_formats"])))
            if row.get("learn_more_url"):
                say("   Learn more: %s" % row["learn_more_url"])
            tail = row.get("source") or row.get("status") or ""
            if tail:
                say("   Source/status: %s" % tail)
    if dropped:
        say("")
        say("DROPPED because their domain is outside your construct:")
        for row in dropped:
            say("  - %s [%s]" % (row["assay_name"], row["entry_id"]))
    if out_of_scope:
        say("")
        say("NOT SHOWN (cell-based; hidden because of --cell-free-only):")
        for row in out_of_scope:
            say("  - %s [%s]" % (row["assay_name"], row["entry_id"]))
    say("")
    say("Ratings use five levels: very low, low, moderate, high, very high.")
    say("Priority scale as I read the reviewer's table (please confirm): P1 develop first; P2 develop next "
        "as follow-up or orthogonal; P3 mechanism or selectivity follow-up; P4 specialized.")
    say("Notes: cell-based assays are labeled; many target classes (GPCRs, ion channels) are screened mainly in cells. Plate format, instruments and reagent availability are not yet applied.")

    name = "m4_%s_%d-%d.txt" % (info["gene"] or accession, start, end)
    out = os.path.join(HERE, name)
    try:
        with open(out, "w", encoding="utf-8") as handle:
            handle.write("\n".join(LOG) + "\n")
        print("\nSaved as: %s" % out)
    except OSError as error:
        print("Could not save (%s)." % error)


if __name__ == "__main__":
    main()
