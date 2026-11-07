"""
M1 data fetch for the Assay Strategy Planner.

What it does: downloads public data for WRN and an EGFR sample and saves it as
files you can open in Excel, plus a short report.

  UniProt  -> the protein's domains, regions and sites, plus database cross-references
  PubChem  -> assays linked to the gene (by three different lookups, compared)
  ChEMBL   -> assays for the protein target, with their descriptions

It only READS public data. Nothing is changed anywhere except new files in a
folder called m1_data next to this script.

How to run (the same way as the M0 test):
    Windows:   python m1_fetch_data.py
    Mac/Linux: python3 m1_fetch_data.py

Optional: add the word wrn to do WRN only (faster, a good first run):
    python m1_fetch_data.py wrn

Downloads are remembered. If you run it again it reuses what it already has.
To start fresh, delete the m1_data folder.

Uses only Python's built-in libraries, so nothing needs installing.
"""

import collections
import csv
import json
import os
import re
import socket
import ssl
import sys
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "m1_data")

# desc_limit: how many PubChem assays get their FULL description downloaded.
# None means all of them. For a big target we take an evenly spaced sample.
TARGETS = [
    {"gene": "WRN", "uniprot": "Q14191", "geneid": 7486, "desc_limit": None},
    {"gene": "EGFR", "uniprot": "P00533", "geneid": 1956, "desc_limit": 40},
]

PAUSE = 0.4          # seconds between calls to the same service (be polite)
RETRY_PAUSE = 3.0    # seconds to wait before retrying a failed call
SUMMARY_BATCH = 25   # how many PubChem assays to ask about at once

HEADERS = {
    "User-Agent": "assay-planner-m1-fetch/0.1 (personal side project)",
    "Accept": "application/json",
}

# UniProt feature types worth keeping (others, such as modified residues and
# secondary structure, are counted in the report but not saved to the CSV).
KEEP_FEATURES = {
    "Domain", "Region", "Repeat", "Motif", "Zinc finger", "DNA binding",
    "Nucleotide binding", "Active site", "Binding site", "Site",
    "Topological domain", "Transmembrane", "Intramembrane", "Signal", "Chain",
    "Compositional bias", "Coiled coil",
}

# Cross-reference databases worth keeping.
KEEP_XREFS = {"Pfam", "InterPro", "SMART", "PROSITE", "PDB", "GO"}

# Crude keyword counts, only to see which technologies show up in the text.
# This is NOT the labeling step (that is M2).
KEYWORDS = collections.OrderedDict([
    ("ADP-Glo", r"adp[- ]?glo"),
    ("Kinase-Glo", r"kinase[- ]?glo"),
    ("Transcreener", r"transcreener"),
    ("HTRF / TR-FRET", r"htrf|tr-?fret|time[- ]resolved (fluorescence|fret)"),
    ("LanthaScreen", r"lanthascreen"),
    ("AlphaScreen / AlphaLISA", r"alpha ?screen|alphalisa"),
    ("Fluorescence polarization", r"fluorescence polari[sz]ation"),
    ("Fluorogenic / quenched probe", r"tamra|quench|fluorogenic|dequench"),
    ("Gel-based / EMSA", r"gel[- ]based|electrophoresis|gel shift|emsa"),
    ("Radiometric", r"33p|32p|radiometric|radioactiv"),
    ("Mobility shift", r"mobility shift|caliper|labchip"),
    ("ELISA", r"elisa"),
    ("SPR", r"surface plasmon|\bspr\b"),
    ("Thermal shift / DSF", r"thermal shift|differential scanning fluorimetry|\bdsf\b"),
    ("ITC", r"isothermal titration"),
    ("Luciferase / luminescence", r"luciferase|luminescen"),
    ("siRNA / RNAi", r"sirna|rnai|shrna"),
    ("Cell-based wording", r"\bcells?\b|cell line|cellular"),
])

LOG = []


def say(text=""):
    print(text)
    LOG.append(text)


# ------------------------------------------------------------ web helpers


def http_get_text(url, timeout=60):
    """The one place that touches the network."""
    request = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8")


def get_json(url, tries=3):
    """Fetch JSON, retrying a few times if the server is busy."""
    last_error = None
    for attempt in range(1, tries + 1):
        try:
            return json.loads(http_get_text(url))
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504):
                raise
            last_error = error
        except (urllib.error.URLError, socket.timeout, ssl.SSLError) as error:
            last_error = error
        if attempt < tries:
            time.sleep(RETRY_PAUSE * attempt)
    raise last_error


def explain(error):
    if isinstance(error, urllib.error.HTTPError):
        return "server error code %s (%s)" % (error.code, error.reason)
    if isinstance(error, urllib.error.URLError):
        return "could not connect (%s)" % error.reason
    return "%s: %s" % (type(error).__name__, error)


def save_json(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle)


def cached_json(url, path):
    """Return (data, downloaded_now). Reuses the file if we already have it."""
    if os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            return json.load(handle), False
    data = get_json(url)
    save_json(path, data)
    return data, True


def write_csv(path, rows, columns):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # utf-8-sig so that Excel shows special characters correctly
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def to_text(value):
    """Turn a string, list or None into one tidy line of text."""
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        value = " ".join(to_text(v) for v in value)
    return " ".join(str(value).split())


def clip(text, limit):
    text = to_text(text)
    return text if len(text) <= limit else text[: limit - 3] + "..."


# ---------------------------------------------------------------- UniProt


def do_uniprot(target):
    acc = target["uniprot"]
    gene = target["gene"]
    url = "https://rest.uniprot.org/uniprotkb/%s.json" % acc
    entry, _ = cached_json(url, os.path.join(DATA, "uniprot_%s.json" % acc))
    features = entry.get("features", [])
    counts = collections.Counter(f.get("type", "?") for f in features)

    rows = []
    for feature in features:
        ftype = feature.get("type", "?")
        if ftype not in KEEP_FEATURES:
            continue
        location = feature.get("location", {})
        start = (location.get("start") or {}).get("value", "")
        end = (location.get("end") or {}).get("value", "")
        description = feature.get("description") or (feature.get("ligand") or {}).get("name", "")
        rows.append({
            "gene": gene, "type": ftype, "start": start, "end": end,
            "description": to_text(description),
        })
    rows.sort(key=lambda r: r["start"] if isinstance(r["start"], int) else 10 ** 9)
    write_csv(os.path.join(DATA, "uniprot_%s_features.csv" % acc), rows,
              ["gene", "type", "start", "end", "description"])

    xrefs = []
    for ref in entry.get("uniProtKBCrossReferences", []):
        database = ref.get("database")
        if database not in KEEP_XREFS:
            continue
        props = {p.get("key"): p.get("value") for p in ref.get("properties", [])}
        if database == "GO":
            term = props.get("GoTerm", "")
            if not term.startswith("F:"):
                continue
            detail = term
        elif database == "PDB":
            detail = "; ".join("%s=%s" % (k, v) for k, v in props.items())
        else:
            detail = props.get("EntryName") or "; ".join(
                "%s=%s" % (k, v) for k, v in props.items())
        xrefs.append({"gene": gene, "database": database, "id": ref.get("id", ""),
                      "detail": to_text(detail)})
    write_csv(os.path.join(DATA, "uniprot_%s_xrefs.csv" % acc), xrefs,
              ["gene", "database", "id", "detail"])

    length = (entry.get("sequence") or {}).get("length", "?")
    say("  UniProt %s: length %s residues, %d annotations (%d saved to CSV)."
        % (acc, length, len(features), len(rows)))
    say("    Annotation types: " + ", ".join(
        "%s %d" % (k, v) for k, v in counts.most_common(12)))
    by_db = collections.Counter(x["database"] for x in xrefs)
    say("    Cross-references saved: " + ", ".join(
        "%s %d" % (k, v) for k, v in sorted(by_db.items())))
    for row in rows:
        if row["type"] in ("Domain", "Repeat", "Topological domain", "Transmembrane"):
            say("    %-18s %5s-%-5s %s" % (row["type"], row["start"], row["end"],
                                            clip(row["description"], 60)))


# ----------------------------------------------------------------- PubChem


def pubchem_aids(kind, value):
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/%s/%s/aids/JSON"
           % (kind, urllib.parse.quote(str(value))))
    try:
        data = get_json(url)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return []
        raise
    return data.get("IdentifierList", {}).get("AID", [])


def pubchem_summaries(gene, aids):
    path = os.path.join(DATA, "pubchem_%s_summaries.json" % gene)
    if os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            return {int(k): v for k, v in json.load(handle).items()}, []
    items, problems = {}, []
    total = len(aids)
    for i in range(0, total, SUMMARY_BATCH):
        batch = aids[i:i + SUMMARY_BATCH]
        url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/summary/JSON"
               % ",".join(str(a) for a in batch))
        try:
            data = get_json(url)
            for item in data["AssaySummaries"]["AssaySummary"]:
                items[int(item["AID"])] = item
        except Exception as error:  # keep going, report at the end
            problems.append("summary batch starting at %d: %s" % (i, explain(error)))
        print("    summaries: %d of %d" % (min(i + SUMMARY_BATCH, total), total))
        time.sleep(PAUSE)
    save_json(path, {str(k): v for k, v in items.items()})
    return items, problems


def parse_description(data):
    """Pull readable text out of PubChem's full assay description."""
    try:
        descr = data["PC_AssayContainer"][0]["assay"]["descr"]
    except (KeyError, IndexError, TypeError):
        return {}
    targets = []
    for t in descr.get("target") or []:
        if isinstance(t, dict) and t.get("name"):
            targets.append(t["name"])
    return {
        "name": to_text(descr.get("name")),
        "description": to_text(descr.get("description")),
        "protocol": to_text(descr.get("protocol")),
        "comment": to_text(descr.get("comment")),
        "targets": "; ".join(targets),
    }


def pick_sample(aids, limit):
    if limit is None or len(aids) <= limit:
        return list(aids)
    step = len(aids) / float(limit)
    return [aids[int(i * step)] for i in range(limit)]


def do_pubchem(target):
    gene = target["gene"]
    methods = collections.OrderedDict()
    lookups = [
        ("genesymbol", gene),
        ("accession", target["uniprot"]),
        ("geneid", target["geneid"]),
    ]
    for kind, value in lookups:
        try:
            methods[kind] = sorted(set(pubchem_aids(kind, value)))
            say("    lookup by %-10s (%s): %d assays" % (kind, value, len(methods[kind])))
        except Exception as error:
            methods[kind] = None
            say("    lookup by %-10s (%s): FAILED, %s" % (kind, value, explain(error)))
        time.sleep(PAUSE)

    all_aids = sorted(set().union(*[set(v) for v in methods.values() if v]))
    if not all_aids:
        say("    No PubChem assays found by any lookup.")
        return {"aids": [], "texts": [], "rows": []}
    if methods.get("genesymbol") is not None and methods.get("accession") is not None:
        a, b = set(methods["genesymbol"]), set(methods["accession"])
        say("    gene symbol vs UniProt accession: %d shared, %d only by symbol, %d only by accession"
            % (len(a & b), len(a - b), len(b - a)))
    say("    %d distinct assays in total; downloading titles..." % len(all_aids))

    summaries, problems = pubchem_summaries(gene, all_aids)
    if summaries:
        keys = sorted(set().union(*[set(s.keys()) for s in summaries.values()]))
        say("    Fields PubChem gives in each summary: " + ", ".join(keys))

    chosen = pick_sample(all_aids, target["desc_limit"])
    say("    downloading full descriptions for %d assays..." % len(chosen))
    full = {}
    for n, aid in enumerate(chosen, 1):
        url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%d/description/JSON" % aid)
        path = os.path.join(DATA, "pubchem_raw", "AID_%d.json" % aid)
        try:
            data, downloaded = cached_json(url, path)
            full[aid] = parse_description(data)
            if downloaded:
                time.sleep(PAUSE)
        except Exception as error:
            problems.append("description for AID %d: %s" % (aid, explain(error)))
        if n % 10 == 0 or n == len(chosen):
            print("    descriptions: %d of %d" % (n, len(chosen)))

    rows, texts = [], []
    for aid in all_aids:
        s = summaries.get(aid, {})
        f = full.get(aid, {})
        found_by = "|".join(k for k, v in methods.items() if v and aid in v)
        title = to_text(s.get("Name")) or f.get("name", "")
        rows.append({
            "aid": aid,
            "source": to_text(s.get("SourceName")),
            "title": title,
            "found_by": found_by,
            "summary_text": clip(s.get("Description"), 1500),
            "full_description": clip(f.get("description", ""), 1500),
            "protocol": clip(f.get("protocol", ""), 1500),
            "targets": f.get("targets", ""),
        })
        texts.append(" ".join([title, to_text(s.get("Description")),
                               f.get("description", ""), f.get("protocol", "")]))
    write_csv(os.path.join(DATA, "pubchem_%s_assays.csv" % gene), rows,
              ["aid", "source", "title", "found_by", "summary_text",
               "full_description", "protocol", "targets"])
    say("    saved %d PubChem rows (%d with full descriptions)." % (len(rows), len(full)))
    for problem in problems[:5]:
        say("    problem: " + problem)
    if len(problems) > 5:
        say("    ... and %d more problems" % (len(problems) - 5))
    if len(rows) <= 40:
        say("    All titles:")
        for row in rows:
            say("      AID %-8s %-28s %s" % (row["aid"], clip(row["source"], 28),
                                              clip(row["title"], 80)))
    return {"aids": all_aids, "texts": texts, "rows": rows}


# ------------------------------------------------------------------ ChEMBL


def chembl_target(target):
    url = ("https://www.ebi.ac.uk/chembl/api/data/target/search.json?q=%s&limit=20"
           % urllib.parse.quote(target["gene"]))
    for t in get_json(url).get("targets", []):
        if t.get("organism") != "Homo sapiens" or t.get("target_type") != "SINGLE PROTEIN":
            continue
        accessions = [c.get("accession") for c in t.get("target_components", [])]
        if target["uniprot"] in accessions:
            return t["target_chembl_id"], t.get("pref_name", "")
    return None, ""


def do_chembl(target, pubchem_aids_found):
    gene = target["gene"]
    target_id, name = chembl_target(target)
    if not target_id:
        say("    ChEMBL: no human single-protein target found for %s." % target["uniprot"])
        return {"texts": []}
    say("    ChEMBL target: %s (%s)" % (name, target_id))
    path = os.path.join(DATA, "chembl_%s_assays_raw.json" % gene)
    if os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            assays = json.load(handle)
    else:
        assays = []
        url = ("https://www.ebi.ac.uk/chembl/api/data/assay.json?target_chembl_id=%s&limit=1000"
               % target_id)
        while url:
            page = get_json(url)
            assays.extend(page.get("assays", []))
            nxt = (page.get("page_meta") or {}).get("next")
            url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
            print("    ChEMBL assays downloaded: %d of %s"
                  % (len(assays), (page.get("page_meta") or {}).get("total_count", "?")))
            time.sleep(PAUSE)
        save_json(path, assays)

    rows, texts = [], []
    for a in assays:
        variant = a.get("variant_sequence") or {}
        rows.append({
            "assay_chembl_id": a.get("assay_chembl_id", ""),
            "assay_type": a.get("assay_type", ""),
            "description": clip(a.get("description"), 1500),
            "bao_label": a.get("bao_label", ""),
            "assay_category": a.get("assay_category", ""),
            "cell_type": a.get("assay_cell_type", ""),
            "organism": a.get("assay_organism", ""),
            "mutation": variant.get("mutation", "") if isinstance(variant, dict) else "",
            "confidence_score": a.get("confidence_score", ""),
            "document": a.get("document_chembl_id", ""),
            "src_id": a.get("src_id", ""),
            "src_assay_id": a.get("src_assay_id", ""),
        })
        texts.append(to_text(a.get("description")))
    write_csv(os.path.join(DATA, "chembl_%s_assays.csv" % gene), rows,
              ["assay_chembl_id", "assay_type", "description", "bao_label",
               "assay_category", "cell_type", "organism", "mutation",
               "confidence_score", "document", "src_id", "src_assay_id"])

    say("    saved %d ChEMBL rows." % len(rows))
    types = collections.Counter(r["assay_type"] for r in rows)
    say("    ChEMBL assay types: " + ", ".join("%s %d" % kv for kv in types.most_common()))
    bao = collections.Counter(r["bao_label"] or "(none)" for r in rows)
    say("    ChEMBL BAO formats: " + ", ".join("%s %d" % kv for kv in bao.most_common(6)))
    mutants = sum(1 for r in rows if r["mutation"])
    say("    Assays on a mutant form of the protein: %d" % mutants)
    pub_like = [r for r in rows if r["description"].lower().startswith("pubchem bioassay")]
    shared = {str(r["src_assay_id"]) for r in pub_like} & {str(a) for a in pubchem_aids_found}
    say("    Assays that are copies of PubChem deposits: %d; of those, %d match an AID "
        "we found in PubChem" % (len(pub_like), len(shared)))
    if len(rows) <= 40:
        say("    All descriptions:")
        for r in rows:
            say("      %s | %s | %s" % (r["assay_chembl_id"], r["assay_type"],
                                         clip(r["description"], 90)))
    return {"texts": texts}


# --------------------------------------------------------------- keywords


def keyword_table(name, texts):
    if not texts:
        return
    say("    Keyword mentions in %s (%d texts):" % (name, len(texts)))
    for label, pattern in KEYWORDS.items():
        regex = re.compile(pattern, re.IGNORECASE)
        hits = sum(1 for t in texts if regex.search(t))
        if hits:
            say("      %-30s %d" % (label, hits))


# -------------------------------------------------------------------- main


def stage(label, function, *args):
    try:
        return function(*args)
    except Exception as error:
        say("  %s FAILED: %s" % (label, explain(error)))
        LOG.append(traceback.format_exc())
        return None


def main():
    say("M1 data fetch for the Assay Strategy Planner")
    say("Run at: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    only_wrn = any(a.lower() == "wrn" for a in sys.argv[1:])
    targets = [t for t in TARGETS if t["gene"] == "WRN"] if only_wrn else TARGETS
    os.makedirs(DATA, exist_ok=True)

    for target in targets:
        say("")
        say("=== %s (UniProt %s) ===" % (target["gene"], target["uniprot"]))
        stage("UniProt", do_uniprot, target)
        say("  PubChem BioAssay:")
        pubchem = stage("PubChem", do_pubchem, target) or {"aids": [], "texts": []}
        say("  ChEMBL:")
        chembl = stage("ChEMBL", do_chembl, target, pubchem["aids"]) or {"texts": []}
        say("  Technology keywords (rough, text only):")
        keyword_table("PubChem text", pubchem["texts"])
        keyword_table("ChEMBL descriptions", chembl["texts"])

    say("")
    say("DONE. Files are in: %s" % DATA)
    report = os.path.join(HERE, "m1_report.txt")
    try:
        with open(report, "w", encoding="utf-8") as handle:
            handle.write("\n".join(LOG) + "\n")
        say("Report saved as: %s" % report)
    except OSError as error:
        say("Could not save the report (%s). Please copy the text above." % error)


if __name__ == "__main__":
    main()
