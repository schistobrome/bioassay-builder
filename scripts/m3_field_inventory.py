"""
M3 field inventory for the Assay Strategy Planner.

Answers one question: WHAT INFORMATION can we actually get about an assay, and
from where? It lists every field the sources give us, how often each is filled,
and which fields hold free text a language model could read.

PART 1 (no internet needed): looks at what M1 already saved
    - every ChEMBL assay record (all fields, how often filled)
    - every PubChem summary and description, split by who deposited it
      (a ChEMBL copy, a BindingDB patent, or another depositor)
    - how ChEMBL assays link to PubChem twins, and whether the twin has more text
    - constructs that crystal structures used (from UniProt's PDB list)
PART 2 (internet, small samples, every answer cached on disk):
    - ChEMBL: activity records, the paper record (document), the target record
    - PubChem: extra endpoints on 3 sample assays (cross-references, targets)
    - Europe PMC: do the papers behind the ChEMBL assays have open full text?
    - RCSB PDB: what a structure entry tells us about the construct
PART 3: checks your GPU (for choosing a bigger language model).

How to run (with the Python inside .venv, in the AssayBuilder folder):
    .\\.venv\\Scripts\\python.exe m3_field_inventory.py
Options:  --offline  (Part 1 and 3 only)     wrn   (WRN only)

Files it creates:
    m3_field_report.txt        paste this back into the chat
    m3_data/field_inventory.csv   every field of every source (open in Excel)
    m3_data/raw/               cached internet answers (reruns call nothing)

It only READS public data. Uses only Python's built-in libraries.
"""

import argparse
import collections
import csv
import glob
import hashlib
import json
import os
import platform
import re
import statistics
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
GENES = {"WRN": "Q14191", "EGFR": "P00533"}
HEADERS = {"User-Agent": "assay-planner-m3/0.1 (personal side project)"}
PAUSE = 0.35
TEXT_MIN_LEN = 60     # a field counts as "free text" if its median length is at least this
LOG = []
INVENTORIES = []


def say(text=""):
    print(text)
    LOG.append(text)


# ---------------------------------------------------------------- helpers


def is_empty(value):
    return value is None or value == "" or value == [] or value == {}


def flatten(obj, path="", out=None, depth=0):
    """List (field path, value) pairs. Lists of plain values are joined into one text."""
    if out is None:
        out = []
    if depth > 7:
        out.append((path, str(obj)[:50]))
    elif isinstance(obj, dict):
        for key, value in obj.items():
            flatten(value, (path + "." + str(key)) if path else str(key), out, depth + 1)
    elif isinstance(obj, list):
        if all(not isinstance(x, (dict, list)) for x in obj):
            out.append((path, " ".join(str(x) for x in obj if x not in (None, ""))))
        else:
            if not obj:
                out.append((path + "[]", ""))
            for item in obj[:50]:
                flatten(item, path + "[]", out, depth + 1)
    else:
        out.append((path, obj))
    return out


class Inventory:
    """Counts how often each field is filled across a set of records."""

    def __init__(self, name):
        self.name = name
        self.n = 0
        self.filled = collections.Counter()
        self.lengths = collections.defaultdict(list)
        self.example = {}
        INVENTORIES.append(self)

    def add(self, record):
        self.n += 1
        seen = set()
        for path, value in flatten(record):
            if is_empty(value):
                self.filled.setdefault(path, 0)
                continue
            text = str(value)
            if path not in seen:
                self.filled[path] += 1
                seen.add(path)
            if len(self.lengths[path]) < 3000:
                self.lengths[path].append(len(text))
            self.example.setdefault(path, " ".join(text.split())[:90])

    def median(self, path):
        lengths = self.lengths.get(path)
        return int(statistics.median(lengths)) if lengths else 0

    def text_fields(self):
        rows = []
        for path, count in self.filled.items():
            if count and self.median(path) >= TEXT_MIN_LEN:
                example = self.example.get(path, "")
                if re.fullmatch(r"[\d\s.,:;/-]+", example):
                    continue
                rows.append((path, 100.0 * count / max(1, self.n), self.median(path)))
        return sorted(rows, key=lambda r: -r[1])

    def show(self, limit=30):
        say("")
        say("  %s: %d records" % (self.name, self.n))
        texty = self.text_fields()
        say("    Free-text fields (what a language model could read):")
        if not texty:
            say("      (none)")
        for path, pct, med in texty[:10]:
            say("      %-52s filled %3d%%, median %5d chars" % (path[-52:], pct, med))
        if limit:
            texty_paths = {p for p, _a, _b in texty}
            others = [(p, c) for p, c in self.filled.items() if p not in texty_paths]
            others.sort(key=lambda r: (-r[1], r[0]))
            say("    Other fields (filled % | example); the CSV has them all:")
            for path, count in others[:limit]:
                say("      %-40s %3d%% | %s" % (path[-40:], 100.0 * count / max(1, self.n),
                                                  self.example.get(path, "")[:44]))
            if len(others) > limit:
                say("      ... and %d more fields" % (len(others) - limit))

    def csv_rows(self):
        for path, count in sorted(self.filled.items()):
            yield {"inventory": self.name, "field": path, "records": self.n, "filled": count,
                   "filled_pct": round(100.0 * count / max(1, self.n), 1),
                   "median_length": self.median(path), "example": self.example.get(path, "")}


def load_json(path):
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def one_line(value):
    if isinstance(value, list):
        value = " ".join(str(x) for x in value)
    return " ".join(str(value or "").split())


# ----------------------------------------------------------------- network


def http_get_text(url, accept="application/json", timeout=60):
    headers = dict(HEADERS)
    headers["Accept"] = accept
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def cached_get(url, raw_dir, accept="application/json"):
    """Fetch a URL once; later runs read it from disk."""
    path = os.path.join(raw_dir, hashlib.sha1(url.encode("utf-8")).hexdigest() + ".txt")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            return handle.read()
    last = None
    for attempt in range(3):
        try:
            text = http_get_text(url, accept)
            os.makedirs(raw_dir, exist_ok=True)
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(text)
            time.sleep(PAUSE)
            return text
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504):
                raise
            last = error
        except (urllib.error.URLError, TimeoutError) as error:
            last = error
        time.sleep(3 * (attempt + 1))
    raise last


def cached_json(url, raw_dir):
    return json.loads(cached_get(url, raw_dir))


def explain(error):
    if isinstance(error, urllib.error.HTTPError):
        return "server answered %s (%s)" % (error.code, error.reason)
    if isinstance(error, urllib.error.URLError):
        return "could not connect (%s)" % error.reason
    return "%s: %s" % (type(error).__name__, error)


# ------------------------------------------------------------------ PART 1


def category(summary):
    source = one_line(summary.get("SourceName")).lower()
    if source == "chembl":
        return "ChEMBL copy"
    if source == "bindingdb":
        return "BindingDB patent"
    return "other depositor"


def part1_gene(gene, m1):
    say("")
    say("=== %s ===" % gene)
    chembl = load_json(os.path.join(m1, "chembl_%s_assays_raw.json" % gene))
    summaries = load_json(os.path.join(m1, "pubchem_%s_summaries.json" % gene)) or {}
    if chembl is None and not summaries:
        say("  No M1 files for %s. Run M1 first." % gene)
        return None, None

    if chembl:
        inv = Inventory("ChEMBL assay record (%s)" % gene)
        for assay in chembl:
            inv.add(assay)
        inv.show(limit=45)

    if summaries:
        by_cat = collections.defaultdict(list)
        for aid, summary in summaries.items():
            by_cat[category(summary)].append(aid)
        say("")
        say("  PubChem assays found for %s, by depositor: %s" % (gene, ", ".join(
            "%s %d" % (k, len(v)) for k, v in sorted(by_cat.items()))))
        for cat_name, aids in sorted(by_cat.items()):
            inv = Inventory("PubChem summary, %s (%s)" % (cat_name, gene))
            for aid in aids:
                inv.add(summaries[aid])
            inv.show(limit=0)
        for cat_name, aids in sorted(by_cat.items()):
            records = []
            for aid in aids:
                data = load_json(os.path.join(m1, "pubchem_raw", "AID_%s.json" % aid))
                if data:
                    records.append(data)
            if records:
                inv = Inventory("PubChem full description JSON, %s (%s)" % (cat_name, gene))
                for record in records:
                    inv.add(record)
                inv.show(limit=25)
            else:
                say("")
                say("  PubChem full description JSON, %s (%s): none saved by M1 (only a sample is saved)"
                    % (cat_name, gene))
    return chembl, summaries


def part1_links(gene, chembl, summaries, m1):
    say("")
    say("  --- Do ChEMBL assays have a PubChem twin with more text? (%s) ---" % gene)
    if not chembl or not summaries:
        say("  Missing data, skipped.")
        return
    src = collections.Counter(str(a.get("src_id")) for a in chembl)
    say("  ChEMBL assays by source id (1 = literature, 7 = PubChem deposit): %s"
        % ", ".join("%s: %d" % kv for kv in sorted(src.items())))
    copies = collections.defaultdict(list)
    for aid, summary in summaries.items():
        if category(summary) == "ChEMBL copy":
            copies[one_line(summary.get("Name"))].append(aid)

    literature = [a for a in chembl if str(a.get("src_id")) == "1"]
    unique = ambiguous = none = 0
    abstract_like = 0
    lengths_pc, lengths_ch, shown = [], [], 0
    for assay in literature:
        title = one_line(assay.get("description"))
        aids = copies.get(title, [])
        if not aids:
            none += 1
            continue
        if len(aids) == 1:
            unique += 1
        else:
            ambiguous += 1
        text = one_line(summaries[aids[0]].get("Description"))
        lengths_pc.append(len(text))
        lengths_ch.append(len(title))
        if text.startswith("Title:"):
            abstract_like += 1
        if shown < 2:
            shown += 1
            say("    example: %s  <->  PubChem AID %s" % (assay.get("assay_chembl_id"), aids[0]))
            say("      ChEMBL text : %s" % title[:110])
            say("      PubChem text: %s" % text[:160])
    say("  Literature assays: %d. PubChem twin found by identical title: %d unique, %d ambiguous "
        "(several twins share a title), %d none." % (len(literature), unique, ambiguous, none))
    if lengths_pc:
        say("  Twin text starts with 'Title: ... Abstract:' (the paper's abstract): %d of %d"
            % (abstract_like, len(lengths_pc)))
        say("  Median text length: ChEMBL %d chars vs PubChem twin %d chars"
            % (statistics.median(lengths_ch), statistics.median(lengths_pc)))

    deposits = [a for a in chembl if str(a.get("src_id")) == "7"]
    have = sum(1 for a in deposits if os.path.exists(os.path.join(
        m1, "pubchem_raw", "AID_%s.json" % one_line(a.get("src_assay_id")))))
    say("  ChEMBL assays that are PubChem deposits: %d; original description already saved for %d."
        % (len(deposits), have))


def pdb_constructs(gene, accession, m1):
    say("")
    say("  --- Constructs that crystal structures used (UniProt -> PDB list) (%s) ---" % gene)
    entry = load_json(os.path.join(m1, "uniprot_%s.json" % accession))
    if not entry:
        say("  No UniProt file saved.")
        return
    domains = []
    for feature in entry.get("features", []):
        if feature.get("type") == "Domain":
            loc = feature.get("location", {})
            domains.append((feature.get("description", "?"), loc.get("start", {}).get("value"),
                            loc.get("end", {}).get("value")))
    ranges = collections.defaultdict(list)
    total = 0
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") != "PDB":
            continue
        total += 1
        props = {p.get("key"): p.get("value") for p in ref.get("properties", [])}
        for match in re.finditer(r"=(\d+)-(\d+)", props.get("Chains", "")):
            ranges[(int(match.group(1)), int(match.group(2)))].append(ref.get("id"))
    say("  %d PDB entries; %d distinct residue ranges." % (total, len(ranges)))
    for (start, end), ids in sorted(ranges.items(), key=lambda kv: -len(kv[1]))[:8]:
        say("    %5d-%-5d  %2d entries (e.g. %s)" % (start, end, len(ids), ", ".join(ids[:3])))
    for name, d_start, d_end in domains:
        if not d_start or not d_end:
            continue
        covering = sum(len(ids) for (s, e), ids in ranges.items() if s <= d_start and e >= d_end)
        say("    domain %-26s %4d-%-4d fully inside %d structure ranges" % (name[:26], d_start, d_end, covering))


# ------------------------------------------------------------------ PART 2

KEYWORDS = collections.OrderedDict([
    ("ADP-Glo", r"adp[- ]?glo"), ("Kinase-Glo", r"kinase[- ]?glo"), ("HTRF/TR-FRET", r"htrf|tr-?fret"),
    ("fluorescence polarization", r"fluorescence polari[sz]ation"), ("gel/EMSA", r"gel[- ]based|polyacrylamide|emsa"),
    ("radiometric", r"\b3[23]p\b"), ("SPR", r"surface plasmon|\bspr\b"), ("thermal shift", r"thermal shift|\bdsf\b"),
    ("ITC", r"isothermal titration"), ("FRET/quencher", r"quencher|tamra|\bfret\b"),
])
RESIDUES = re.compile(r"\b\d{1,4}\s*(?:to|-|–)\s*\d{1,4}\s*(?:residues|aa|amino acids)\b|\b(?:residues|aa)\s*\d{1,4}\s*(?:to|-|–)\s*\d{1,4}\b", re.I)


def probe(label, url, raw_dir, inv, pick=None):
    """Fetch one JSON answer, add it to an inventory, and say what happened."""
    try:
        data = cached_json(url, raw_dir)
    except Exception as error:
        say("    %-34s FAILED: %s" % (label, explain(error)))
        return None
    records = pick(data) if pick else [data]
    for record in records:
        inv.add(record)
    say("    %-34s ok (%d record%s)" % (label, len(records), "" if len(records) == 1 else "s"))
    return data


def part2(genes, m1, out_dir):
    raw = os.path.join(out_dir, "raw")
    say("")
    say("=== PART 2: online probes (small samples; answers are cached) ===")
    chembl_base = "https://www.ebi.ac.uk/chembl/api/data/"

    for gene in genes:
        assays = load_json(os.path.join(m1, "chembl_%s_assays_raw.json" % gene)) or []
        if not assays:
            continue
        say("")
        say("  -- ChEMBL: %s --" % gene)
        step = max(1, len(assays) // 10)
        sample = assays[::step][:10]
        inv = Inventory("ChEMBL activity record (%s, 5 per assay for %d assays)" % (gene, len(sample)))
        for assay in sample:
            url = chembl_base + "activity.json?assay_chembl_id=%s&limit=5" % assay.get("assay_chembl_id")
            probe("activities of " + str(assay.get("assay_chembl_id")), url, raw, inv,
                  lambda d: d.get("activities", []))
        inv.show(limit=40)

        target_ids = sorted({a.get("target_chembl_id") for a in assays if a.get("target_chembl_id")})
        if target_ids:
            tinv = Inventory("ChEMBL target record (%s)" % gene)
            probe("target " + target_ids[0], chembl_base + "target/%s.json" % target_ids[0], raw, tinv)
            tinv.show(limit=25)

        doc_ids = sorted({a.get("document_chembl_id") for a in assays if a.get("document_chembl_id")})
        if gene == "EGFR":
            doc_ids = doc_ids[:: max(1, len(doc_ids) // 12)][:12]
        dinv = Inventory("ChEMBL document (the paper) record (%s, %d papers)" % (gene, len(doc_ids)))
        docs = []
        for doc_id in doc_ids:
            data = probe("document " + doc_id, chembl_base + "document/%s.json" % doc_id, raw, dinv)
            if data:
                docs.append(data)
        dinv.show(limit=25)
        say("")
        say("  Papers behind these assays: %d with a PubMed id, %d with a DOI, %d patents or other."
            % (sum(1 for d in docs if d.get("pubmed_id")), sum(1 for d in docs if d.get("doi")),
               sum(1 for d in docs if not d.get("pubmed_id"))))
        europepmc(docs, raw)

    summaries = load_json(os.path.join(m1, "pubchem_%s_summaries.json" % genes[0])) or {}
    pubchem_extra(summaries, m1, raw)
    for gene in genes:
        entry = load_json(os.path.join(m1, "uniprot_%s.json" % GENES[gene]))
        if entry:
            rcsb(entry, raw)
            break


def europepmc(docs, raw):
    say("")
    say("  -- Europe PMC: are the papers' methods readable? --")
    inv = Inventory("Europe PMC article record")
    opened = 0
    for doc in docs[:12]:
        pmid = doc.get("pubmed_id")
        if not pmid:
            continue
        url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%s&format=json&resultType=core&pageSize=1"
               % urllib.parse.quote("EXT_ID:%s AND SRC:MED" % pmid))
        try:
            results = cached_json(url, raw)["resultList"]["result"]
        except Exception as error:
            say("    PMID %-9s FAILED: %s" % (pmid, explain(error)))
            continue
        if not results:
            say("    PMID %-9s not found" % pmid)
            continue
        record = results[0]
        inv.add(record)
        pmcid, free = record.get("pmcid"), record.get("isOpenAccess") == "Y"
        if not (pmcid and free):
            say("    PMID %-9s abstract only (open full text: no)" % pmid)
            continue
        try:
            xml = cached_get("https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML" % pmcid,
                             raw, accept="application/xml")
        except Exception as error:
            say("    PMID %-9s open access but full text FAILED: %s" % (pmid, explain(error)))
            continue
        text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml))
        found = [name for name, pattern in KEYWORDS.items() if re.search(pattern, text, re.I)]
        opened += 1
        say("    PMID %-9s full text %6d chars | residue ranges: %d | methods words: %s"
            % (pmid, len(text), len(RESIDUES.findall(text)), ", ".join(found) or "none"))
    say("  Full text readable for %d of %d papers checked." % (opened, min(12, len(docs))))
    inv.show(limit=15)


def pubchem_extra(summaries, m1, raw):
    say("")
    say("  -- PubChem: extra endpoints on sample assays --")
    picks = {}
    for aid, summary in summaries.items():
        picks.setdefault(category(summary), aid)
    if os.path.exists(os.path.join(m1, "pubchem_raw", "AID_651768.json")):
        picks["NCGC original (via ChEMBL)"] = "651768"
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/%s/JSON"
    endpoints = [("xrefs", "xrefs/PubMedID"), ("xrefs-patent", "xrefs/PatentID"),
                 ("targets", "targets/ProteinName,GeneSymbol,GeneID")]
    for name, aid in picks.items():
        for label, path in endpoints:
            inv = Inventory("PubChem %s of AID %s (%s)" % (label, aid, name))
            data = probe("%s AID %s (%s)" % (label, aid, name[:14]), base % (aid, path), raw, inv)
            if data is not None:
                say("        keys: %s" % ", ".join(p for p, _v in flatten(data)[:6])[:150])


def rcsb(entry, raw):
    say("")
    say("  -- RCSB PDB: what a structure entry says about its construct --")
    ids = []
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "PDB" and len(ids) < 2:
            ids.append(ref.get("id"))
    for pdb_id in ids:
        inv = Inventory("RCSB entry %s" % pdb_id)
        probe("entry " + pdb_id, "https://data.rcsb.org/rest/v1/core/entry/%s" % pdb_id, raw, inv)
        inv.show(limit=0)
        inv2 = Inventory("RCSB polymer entity %s/1" % pdb_id)
        probe("polymer entity %s/1" % pdb_id, "https://data.rcsb.org/rest/v1/core/polymer_entity/%s/1" % pdb_id, raw, inv2)
        inv2.show(limit=18)


# ------------------------------------------------------------------ PART 3


def hardware():
    say("")
    say("=== PART 3: your computer (for choosing a bigger language model) ===")
    say("  %s, Python %s" % (platform.platform(), sys.version.split()[0]))
    try:
        result = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total,memory.used,driver_version",
                                 "--format=csv,noheader"], capture_output=True, text=True, timeout=20)
        say("  nvidia-smi: %s" % (result.stdout.strip() or result.stderr.strip() or "no output"))
    except Exception as error:
        say("  nvidia-smi not available (%s)" % type(error).__name__)
    try:
        import torch
        say("  torch %s; CUDA available: %s" % (torch.__version__, torch.cuda.is_available()))
        if torch.cuda.is_available():
            props = torch.cuda.get_device_properties(0)
            say("  GPU seen by torch: %s, %.1f GB" % (props.name, props.total_memory / 1e9))
        elif "+cpu" in torch.__version__:
            say("  This torch is the CPU-only build, so the GPU is NOT being used. A CUDA build is needed.")
    except Exception as error:
        say("  torch not importable here (%s). Run with the .venv Python for this check." % type(error).__name__)


# -------------------------------------------------------------------- main


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genes", nargs="*")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--m1", default=os.path.join(HERE, "m1_data"))
    parser.add_argument("--out", default=os.path.join(HERE, "m3_data"))
    args = parser.parse_args()
    genes = [g.upper() for g in args.genes if g.upper() in GENES] or list(GENES)
    os.makedirs(args.out, exist_ok=True)

    say("M3 field inventory for the Assay Strategy Planner")
    say("Run at: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    say("")
    say("=== PART 1: what M1 already saved (no internet) ===")
    for gene in genes:
        chembl, summaries = part1_gene(gene, args.m1)
        if chembl is not None or summaries:
            part1_links(gene, chembl, summaries, args.m1)
        pdb_constructs(gene, GENES[gene], args.m1)

    if not args.offline:
        try:
            part2(genes, args.m1, args.out)
        except Exception as error:
            say("PART 2 stopped early: %s" % explain(error))
    hardware()

    with open(os.path.join(args.out, "field_inventory.csv"), "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["inventory", "field", "records", "filled",
                                                     "filled_pct", "median_length", "example"])
        writer.writeheader()
        for inv in INVENTORIES:
            writer.writerows(inv.csv_rows())
    say("")
    say("DONE. Every field of every source: %s" % os.path.join(args.out, "field_inventory.csv"))
    report = os.path.join(HERE, "m3_field_report.txt")
    try:
        with open(report, "w", encoding="utf-8") as handle:
            handle.write("\n".join(LOG) + "\n")
        say("Report saved as: %s (please paste it back into the chat)" % report)
    except OSError as error:
        say("Could not save the report (%s). Please copy the text above." % error)


if __name__ == "__main__":
    main()
