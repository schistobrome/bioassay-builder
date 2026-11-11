#!/usr/bin/env python3
"""
One-hit database build, version 2: every human UniProt protein -> ChEMBL + PubChem -> labels -> usage rows.

    python scripts/build_database.py                      # everything (resumes where it stopped)
    python scripts/build_database.py --limit 50           # first 50 proteins only (good first test)
    python scripts/build_database.py --genes EGFR,WRN,BRD4
    python scripts/build_database.py --max-minutes 120    # stop politely after 2 hours; rerun to continue
    python scripts/build_database.py --merge              # also add new rows to data/usage_stats.csv
    python scripts/build_database.py --status             # progress, including parked proteins

WHAT IS NEW IN VERSION 2: no single protein can sink the run
  * Every protein gets a TIME BUDGET that is enforced inside every network call (not only between proteins).
  * PASS 1 ("light"): every protein gets a short budget (default 6 minutes) and a size cap (default 2500 PubChem
    assays). Normal proteins finish. A protein that is too big or too slow is PARKED, not failed: everything it
    already downloaded is kept, and the script moves on to the next protein.
  * PASS 2 ("heavy"): the parked proteins are then worked on, smallest remaining job first, each with a larger
    budget (default 40 minutes) and cap (default 15000 assays; 0 = no cap). A protein that runs out of time is
    parked again and simply continues next run from where it stopped.
  * --max-minutes is a hard limit for the whole run, also inside a protein.
  * Downloads are saved in small chunks (summaries, ChEMBL pages, labels), so nothing is ever downloaded twice.
  * One bad PubChem assay record is isolated (batches are halved) and skipped instead of failing the protein.
  * A server that keeps failing triggers a cool-down; a protein that fails 3 runs in a row is quarantined.
  * Cache layout is unchanged from version 1 (build_cache/), so your existing runs are reused as they are.

What it does, in order:
  A  list      Downloads the reviewed (canonical) human UniProtKB proteins (about 20,000).
  B  filter    Keeps proteins that have a ChEMBL single-protein target (use --include-no-chembl to try the rest).
  C  per protein
       fetch   UniProt entry, ChEMBL assays, PubChem assays (chunked and budgeted).
       label   Term-table synonyms first, then the small local meaning model for gaps (m2.label_records),
               in slices of 300 records so that a huge protein can be paused and continued.
       summary Compact per-protein summary: classes, campaigns (ChEMBL copies of PubChem deposits are merged
               with the PubChem original so one campaign is counted once), function and technology labels.
  D  aggregate Usage rows per protein and per function class, effect counts, readout suggestions, gaps.
  E  merge     (only with --merge) Adds the rows to data/usage_stats.csv; never touches rows you made by hand.

Updating: --refresh-older-than 90, --refresh EGFR,WRN. Changing term_table.csv makes only labeling and summary
run again (no downloads). Standard library only, plus the optional model that m2_label.py already uses.
"""

import argparse
import collections
import contextlib
import csv
import glob
import hashlib
import io
import json
import os
import re
import shutil
import signal
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timedelta

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import m1_fetch_data as m1      # noqa: E402
import m2_label as m2           # noqa: E402
import m4_recommend as m4       # noqa: E402

CACHE = os.path.join(ROOT, "build_cache")
OUTDIR = os.path.join(ROOT, "build_output")
STATE = os.path.join(CACHE, "state")
SUMMARY = os.path.join(CACHE, "summary")
SUMMARY_BDB = os.path.join(CACHE, "summary_bdb")      # written by bindingdb_label.py (optional)
LABELS = os.path.join(CACHE, "m2_data")
LOGS = os.path.join(CACHE, "logs")
DATA = os.path.join(ROOT, "data")
UA = {"User-Agent": "bioassay-builder/0.2 (personal side project)", "Accept": "application/json"}
LABEL_COLUMNS = ["id", "source", "campaign", "text_level", "in_v1_scope", "context", "role", "effect", "mutant",
                 "function_labels", "function_method", "technology_labels", "technology_method",
                 "construct_start", "construct_end", "construct_note", "chembl_type", "bao", "pubchem_aid",
                 "title", "function_evidence", "technology_evidence", "construct_evidence"]
LABEL_FIELDS = ["context", "role", "effect", "mutant", "in_v1_scope", "construct_start", "construct_end",
                "construct_note", "construct_evidence", "function_labels", "function_method", "function_evidence",
                "technology_labels", "technology_method", "technology_evidence"]
USAGE_COLUMNS = ["stat_id", "entry_id", "scope_type", "scope_id", "stage", "effect", "readout", "n_campaigns_total",
                 "n_campaigns_classifiable", "n_campaigns_using", "n_literature_total", "n_literature_using",
                 "source", "query", "retrieved_on", "method_version", "notes"]
PUBCHEM_COLUMNS = ["aid", "source", "title", "found_by", "summary_text", "full_description", "protocol", "targets"]
CHEMBL_COLUMNS = ["assay_chembl_id", "assay_type", "description", "bao_label", "assay_category", "cell_type",
                  "organism", "mutation", "confidence_score", "document", "src_id", "src_assay_id"]
STAGES = ["fetch", "label", "summary"]
STATE_VERSION = 2
MAX_ATTEMPTS = 3          # failed runs before a protein is quarantined (partial downloads: accepted after 3 tries)
AUTO_PREFIXES = ("bd-", "pipeline-")      # method_version prefixes this tooling may replace in usage_stats.csv
LABEL_SLICE = 300
SUMMARY_BATCH = 25
CHEMBL_PAGE = 200
AUTO_ROW_NOTE = "auto-generated; review before use"


class NetworkDown(Exception):
    pass


class BudgetExceeded(BaseException):
    """Raised inside the network door / loops when time is up. BaseException so no 'except Exception' eats it."""

    def __init__(self, kind, reason):
        super().__init__(reason)
        self.kind, self.reason = kind, reason          # kind: "protein" or "global"


class ParkedCap(Exception):
    def __init__(self, n, cap):
        super().__init__("%d PubChem assays > light cap %d" % (n, cap))
        self.n, self.cap = n, cap


# ---------------------------------------------------------------- small helpers

def now():
    return datetime.now().replace(microsecond=0).isoformat()


def log(text=""):
    print(text, flush=True)
    if log.handle:
        log.handle.write(text + "\n")
        log.handle.flush()


log.handle = None


def split(cell):
    return [x.strip() for x in re.split(r"[;|]", cell or "") if x.strip()]


def safe(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def write_json(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(data, handle)
    os.replace(tmp, path)          # a crash never leaves half a file


def read_json(path, default=None):
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return default


def read_csv(path):
    try:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            return list(csv.DictReader(handle))
    except OSError:
        return []


def write_csv(path, rows, columns):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def append_lines(path, objs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", encoding="utf-8") as handle:
        for obj in objs:
            handle.write(json.dumps(obj) + "\n")


def read_lines(path):
    out = []
    try:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                try:
                    out.append(json.loads(line))
                except ValueError:
                    pass                      # a half-written last line from a crash
    except OSError:
        pass
    return out


def file_hash(path):
    try:
        with open(path, "rb") as handle:
            return hashlib.sha1(handle.read()).hexdigest()[:12]
    except OSError:
        return ""


# ---------------------------------------------------------------- budgets, breaker, the one network door

class Budget:
    """Two clocks: one for the protein being worked on, one for the whole run (--max-minutes)."""

    def __init__(self):
        self.protein_deadline = None
        self.global_deadline = None

    def start_protein(self, seconds):
        self.protein_deadline = time.time() + seconds if seconds else None

    def end_protein(self):
        self.protein_deadline = None

    def remaining(self):
        ds = [d for d in (self.protein_deadline, self.global_deadline) if d]
        return min(ds) - time.time() if ds else 1e9

    def check(self):
        t = time.time()
        if self.global_deadline and t >= self.global_deadline:
            raise BudgetExceeded("global", "run time limit (--max-minutes) reached")
        if self.protein_deadline and t >= self.protein_deadline:
            raise BudgetExceeded("protein", "protein time budget used up")

    def sleep(self, seconds):
        self.check()
        time.sleep(max(0.0, min(seconds, self.remaining())))
        self.check()


BUDGET = Budget()


class Breaker:
    """If a host keeps failing, stop hammering it: cool down (inside the budget), then try again."""

    def __init__(self):
        self.fails, self.opens, self.until = {}, {}, {}

    def before(self, host):
        wait = self.until.get(host, 0) - time.time()
        if wait > 0:
            BUDGET.sleep(wait)

    def ok(self, host):
        self.fails[host] = 0
        self.opens[host] = 0

    def bad(self, host):
        self.fails[host] = self.fails.get(host, 0) + 1
        if self.fails[host] >= 4:
            n = self.opens.get(host, 0)
            self.until[host] = time.time() + min(300, 20 * (2 ** n))
            self.opens[host] = n + 1
            self.fails[host] = 0
            if n + 1 >= 4:
                raise NetworkDown("%s keeps failing" % host)


BREAKER = Breaker()
PAUSE = [0.4]
TIMEOUT = [45.0]         # seconds per request (--timeout)
PACE = [1.0]             # multiplier on --pause: rises when a server says "slow down", relaxes when calls succeed
WORK = [0]               # rough count of calls still to do for the protein in progress (used to order the heavy pass)


def _raw_get(url, timeout):
    """The only function that really touches the internet (tests replace it)."""
    request = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace"), dict(response.headers)


def http_get(url, timeout=None):
    """Return (text, headers). Budget- and breaker-aware."""
    host = urllib.parse.urlsplit(url).netloc
    BUDGET.check()
    BREAKER.before(host)
    limit = max(5.0, min(float(timeout or TIMEOUT[0]), BUDGET.remaining()))
    try:
        result = _raw_get(url, limit)
    except urllib.error.HTTPError as error:
        if error.code in (429, 502, 503, 504):          # 500 is often one bad record, not a sick server
            PACE[0] = min(10.0, PACE[0] * 1.5)          # slow down for everything that follows
            try:
                wait = float((error.headers or {}).get("Retry-After", ""))
            except (TypeError, ValueError):
                wait = 0.0
            if wait:
                BUDGET.sleep(min(wait, 60.0))           # the server told us how long to wait
            BREAKER.bad(host)
        else:
            BREAKER.ok(host)
        raise
    except (urllib.error.URLError, OSError):
        PACE[0] = min(10.0, PACE[0] * 1.3)
        BREAKER.bad(host)
        raise
    BREAKER.ok(host)
    PACE[0] = max(1.0, PACE[0] * 0.97)
    return result


def get_json_h(url, tries=3, pause=3.0):
    last = None
    for attempt in range(1, tries + 1):
        try:
            text, headers = http_get(url)
            return json.loads(text), headers
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504):
                raise
            last = error
        except (urllib.error.URLError, OSError) as error:
            last = error
        if attempt < tries:
            BUDGET.sleep(pause * attempt)
    raise last


def cached_json(url, path, tries=3):
    data = read_json(path)
    if data is not None:
        return data
    data, _ = get_json_h(url, tries=tries)
    write_json(path, data)
    return data


def route_m1_through_http_get():
    """Any leftover m1 call uses the same door."""
    m1.http_get_text = lambda url, timeout=60: http_get(url, timeout)[0]
    m1.HEADERS.update(UA)
    m1.RETRY_PAUSE = 1.0


# one-off bad records: a PubChem assay that always breaks gets skipped after two failed attempts
def bad_path():
    return os.path.join(CACHE, "bad_aids.json")


def bad_count(tag):
    return (read_json(bad_path(), {}) or {}).get(tag, 0)


def bad_add(tag):
    data = read_json(bad_path(), {}) or {}
    data[tag] = data.get(tag, 0) + 1
    write_json(bad_path(), data)


# ---------------------------------------------------------------- extra terms for the wider playbook
# The original term table knew kinases, helicases and nucleases. These terms let the labeler recognise the
# newer playbook classes. They are appended to term_table.csv only if missing (your edits are never overwritten).
# (axis, term, synonyms separated by |, example sentences separated by |, suppresses, applies_to)
EXTRA_TERMS = [
    ("function", "Histone deacetylase activity", "hdac|histone deacetylase|deacetylase activity|sirtuin|sirt1|sirt2|sirt3|hdac1|hdac2|hdac3|hdac6|hdac8", "Inhibition of recombinant HDAC1 deacetylase activity using a fluorogenic acetyl-lysine peptide", "", "any"),
    ("function", "Histone methyltransferase activity", "histone methyltransferase|methyltransferase activity|hmt|ezh2|prmt1|prmt5|dot1l|g9a|smyd2|setd2|nsd2|sam methyl transfer", "Inhibition of EZH2 methyltransferase activity measured by methylation of a histone H3 peptide", "", "any"),
    ("function", "Histone demethylase activity", "demethylase|histone demethylase|lsd1|kdm1a|kdm4|kdm5|jmjd|jumonji", "Inhibition of LSD1 demethylase activity in a peroxidase-coupled assay", "", "any"),
    ("function", "Bromodomain binding", "bromodomain|brd4|brd2|brd3|brdt|bet bromodomain|acetyl lysine peptide displacement|acetylated peptide binding|kac peptide", "Displacement of an acetylated histone H4 peptide from the BRD4 bromodomain", "", "any"),
    ("function", "Protease activity", "protease|proteinase|peptidase|caspase|cathepsin|mmp|metalloproteinase|cleavage of a fluorogenic peptide|protease inhibition", "Inhibition of cathepsin K cleavage of a fluorogenic peptide substrate", "", "any"),
    ("function", "Phosphatase activity", "phosphatase|ptp1b|shp2|dephosphorylation|pp2a|dusp|pten", "Inhibition of PTP1B dephosphorylation of pNPP or a phosphopeptide", "", "any"),
    ("function", "Ubiquitination activity", "ubiquitination|ubiquitin ligase|e3 ligase|e2 enzyme|autoubiquitination|ubiquitin transfer|polyubiquitin", "Inhibition of E3 ligase auto-ubiquitination in a TR-FRET assay", "", "any"),
    ("function", "Deubiquitinase activity", "deubiquitinase|deubiquitinating|dub|usp7|usp1|uchl1|ubiquitin-rhodamine|ubiquitin aminomethylcoumarin", "Inhibition of USP7 cleavage of ubiquitin-rhodamine110", "", "any"),
    ("function", "Topoisomerase activity", "topoisomerase|topo i|topo ii|dna relaxation|supercoiled dna relaxation|decatenation|kdna decatenation|cleavage complex", "Inhibition of topoisomerase II kDNA decatenation", "", "any"),
    ("function", "Protein-protein interaction", "protein-protein interaction|ppi|disruption of the interaction|interaction between|bcl-2 family binding|mdm2 p53", "Disruption of the MDM2-p53 interaction measured by fluorescence polarization", "", "any"),
    ("function", "Nuclear receptor coactivator recruitment", "coactivator recruitment|coactivator peptide|nuclear receptor coactivator|src-1 peptide|ligand binding domain coactivator|lbd coactivator", "Recruitment of a fluorescent coactivator peptide to the ligand-binding domain measured by TR-FRET", "", "any"),
    ("function", "RNA-protein binding", "rna binding|rna-protein|rna protein interaction|rna recognition|riboswitch|rna aptamer", "Inhibition of protein binding to a fluorescent RNA oligonucleotide", "", "any"),
    ("function", "G-protein activation (GTPgammaS)", "gtpgammas|gtp gamma s|35s gtp|eu-gtp|g protein activation|g-protein activation", "Agonist-stimulated [35S]GTPgammaS binding in receptor membranes", "", "any"),
    ("function", "Receptor radioligand binding", "radioligand binding|radioligand displacement|competition binding membranes|ligand binding membranes|tritiated ligand|3h ligand|3h labeled|125i ligand|scintillation proximity binding", "Displacement of [3H]-ligand from membranes expressing the receptor", "", "any"),
    ("function", "cAMP signaling", "camp accumulation|camp assay|cyclic amp|intracellular camp|camp levels|camp production|camp htrf|glosensor camp|forskolin", "Agonist-induced cAMP accumulation in HEK293 cells expressing the receptor", "", "cell-based"),
    ("function", "Calcium mobilization", "calcium mobilization|calcium flux|intracellular calcium|calcium release|fluo-4|fluo4|flipr calcium|fura-2|ca2+ mobilization|calcium signaling assay", "Agonist-induced intracellular calcium mobilization measured with Fluo-4", "", "cell-based"),
    ("function", "Beta-arrestin recruitment", "beta-arrestin|b-arrestin|beta arrestin|arrestin recruitment|pathhunter|tango assay|arrestin translocation", "Beta-arrestin recruitment measured by enzyme fragment complementation", "", "cell-based"),
    ("function", "Receptor internalization", "receptor internalization|receptor internalisation|receptor trafficking|receptor endocytosis|surface receptor loss", "Agonist-induced receptor internalization measured by high-content imaging", "", "cell-based"),
    ("function", "Reporter gene expression", "reporter gene|reporter assay|luciferase reporter|cre-luc|cre luciferase|nfat luciferase|sre luciferase|gal4 luciferase|beta-lactamase reporter|transactivation assay", "CRE-luciferase reporter response in cells expressing the receptor", "", "cell-based"),
    ("function", "Hormone secretion", "insulin secretion|glucose-stimulated insulin|gsis|hormone secretion|glp-1 secretion|secretion assay", "Glucose-stimulated insulin secretion in INS-1 cells", "", "cell-based"),
    ("function", "Ion channel activity", "ion channel|patch clamp|patch-clamp|membrane potential|thallium flux|ion flux|herg|sodium channel|potassium channel|calcium channel|nav1|kv11|channel current|automated electrophysiology|qpatch|patchxpress|synchropatch", "Inhibition of hERG current measured by automated patch clamp", "", "cell-based"),
    ("function", "Transporter activity", "transporter activity|substrate uptake|cellular uptake of|efflux ratio|p-glycoprotein|p-gp|abcb1|bcrp|mrp1|slc|vesicular transport|bidirectional transport|caco-2 transport|calcein-am efflux|transporter atpase", "Inhibition of P-gp mediated efflux in a calcein-AM assay", "", "cell-based"),
    ("function", "Cellular protein degradation", "protein degradation|dc50|dmax|hibit|protein level reduction|targeted protein degradation|protac|molecular glue", "Degradation of the target protein measured by HiBiT lytic detection", "", "cell-based"),
    ("function", "Cellular histone modification", "histone mark|h3k27me3|h3k4me|h3k9|h3k27ac|histone acetylation level|histone methylation level|cellular histone", "Reduction of cellular H3K27me3 measured by AlphaLISA", "", "cell-based"),
    ("function", "Autophagic flux", "autophagic flux|lc3 puncta|mcherry-gfp-lc3|autophagy flux|autophagosome|lc3 turnover", "Autophagic flux measured by tandem mCherry-GFP-LC3 imaging", "", "cell-based"),
    ("technology", "Luminescent HDAC assay", "hdac-glo|hdac glo|luminogenic hdac|luminescent hdac", "HDAC-Glo I/II luminescent assay", "", "any"),
    ("technology", "Fluorescent calcium dye", "fluo-4|fluo4|fura-2|calcium dye|fluorescent calcium|flipr|fdss|calcium 5 dye|calcium 6", "FLIPR calcium dye readout", "", "cell-based"),
    ("technology", "Enzyme fragment complementation", "enzyme fragment complementation|efc|pathhunter|hibit complementation|beta-galactosidase complementation", "PathHunter enzyme fragment complementation readout", "", "any"),
    ("technology", "BRET / NanoBRET", "bret|nanobret|bioluminescence resonance energy transfer|nanoluc", "NanoBRET target engagement", "", "cell-based"),
    ("technology", "Automated patch clamp", "automated patch clamp|patch clamp|qpatch|patchxpress|syncropatch|synchropatch|ionworks|ionflux|whole-cell recording", "Automated patch clamp recording", "", "cell-based"),
    ("technology", "Membrane potential dye", "membrane potential dye|fmp dye|flipr membrane potential|voltage sensitive dye|dibac|fluorescent membrane potential", "FLIPR membrane potential dye", "", "cell-based"),
    ("technology", "Radioligand binding", "radioligand|scintillation proximity|spa beads|filtration binding|tritiated|3h ligand|125i ligand", "Scintillation proximity radioligand binding", "", "any"),
    ("technology", "High-content imaging", "high content|high-content|operetta|incell analyzer|cellomics|immunofluorescence", "High-content imaging readout", "", "cell-based"),
]


def ensure_terms():
    """Append any of EXTRA_TERMS that the term table does not have yet. Returns how many were added."""
    if not os.path.exists(m2.TERM_TABLE):
        m2.write_seed_table()
    rows = read_csv(m2.TERM_TABLE)
    have = {(r["axis"].strip().lower(), r["term"].strip().lower()) for r in rows}
    new = [t for t in EXTRA_TERMS if (t[0], t[1].lower()) not in have]
    if not new:
        return 0
    cols = list(rows[0].keys()) if rows else m2.TERM_COLUMNS
    for axis, term, syn, ex, sup, applies in new:
        rows.append({"axis": axis, "term": term, "synonyms": syn, "examples": ex, "suppresses": sup,
                     "applies_to": applies, "notes": "added by build_database.py"})
    write_csv(m2.TERM_TABLE, rows, cols)
    return len(new)


# ---------------------------------------------------------------- A: protein list, B: ChEMBL targets

def next_link(headers):
    link = headers.get("Link") or headers.get("link") or ""
    match = re.search(r'<([^>]+)>;\s*rel="next"', link)
    return match.group(1) if match else None


def load_protein_list(refresh):
    final = os.path.join(CACHE, "uniprot_human_reviewed.csv")
    if os.path.exists(final) and not refresh:
        return read_csv(final)
    folder = os.path.join(CACHE, "uniprot_list")
    if refresh and os.path.isdir(folder):
        shutil.rmtree(folder)
    url = ("https://rest.uniprot.org/uniprotkb/search?query=%s&format=json&size=500&fields=%s"
           % (urllib.parse.quote("(organism_id:9606) AND (reviewed:true)"),
              "accession,gene_primary,protein_name,length"))
    rows, page = [], 0
    while url:
        page += 1
        path = os.path.join(folder, "page_%04d.json" % page)
        saved = read_json(path)
        if saved is None:
            data, headers = get_json_h(url)
            saved = {"results": data.get("results", []), "next": next_link(headers)}
            write_json(path, saved)
        for item in saved["results"]:
            genes = item.get("genes") or [{}]
            name = ""
            desc = item.get("proteinDescription") or {}
            for key in ("recommendedName", "submissionNames"):
                node = desc.get(key)
                node = node[0] if isinstance(node, list) and node else node
                if node and (node.get("fullName") or {}).get("value"):
                    name = node["fullName"]["value"]
                    break
            rows.append({"accession": item["primaryAccession"],
                         "gene": (genes[0].get("geneName") or {}).get("value", ""),
                         "name": name, "length": (item.get("sequence") or {}).get("length", "")})
        log("  protein list: page %d, %d proteins so far" % (page, len(rows)))
        url = saved["next"]
        time.sleep(0.2)
    write_csv(final, rows, ["accession", "gene", "name", "length"])
    return rows


def load_chembl_targets(refresh):
    final = os.path.join(CACHE, "chembl_human_targets.json")
    if os.path.exists(final) and not refresh:
        return read_json(final, {})
    folder = os.path.join(CACHE, "chembl_targets")
    if refresh and os.path.isdir(folder):
        shutil.rmtree(folder)
    url = ("https://www.ebi.ac.uk/chembl/api/data/target.json?organism=%s&target_type=%s&limit=1000"
           % (urllib.parse.quote("Homo sapiens"), urllib.parse.quote("SINGLE PROTEIN")))
    mapping, page = {}, 0
    while url:
        page += 1
        path = os.path.join(folder, "page_%04d.json" % page)
        data = read_json(path)
        if data is None:
            data, _ = get_json_h(url)
            write_json(path, data)
        for target in data.get("targets", []):
            for comp in target.get("target_components", []):
                if comp.get("accession"):
                    mapping[comp["accession"]] = target["target_chembl_id"]
        nxt = (data.get("page_meta") or {}).get("next")
        url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
        log("  ChEMBL targets: page %d, %d accessions so far" % (page, len(mapping)))
        time.sleep(0.2)
    write_json(final, mapping)
    return mapping


# ---------------------------------------------------------------- per-protein state

def state_path(acc):
    return os.path.join(STATE, safe(acc) + ".json")


def load_state(acc):
    st = read_json(state_path(acc), None)
    if not isinstance(st, dict) or "stages" not in st:
        st = {"acc": acc, "stages": {}, "v": STATE_VERSION}
    if st.get("v") != STATE_VERSION:
        # states written by version 1 (or by the old timeout behaviour): give failed steps a fresh start
        for info in st["stages"].values():
            if info.get("status") == "failed":
                info["attempts"] = 0
        st["v"] = STATE_VERSION
    return st


def save_state(st):
    write_json(state_path(st["acc"]), st)


def stage_done(st, name):
    return (st["stages"].get(name) or {}).get("status") == "done"


def stage_status(st, name):
    return (st["stages"].get(name) or {}).get("status", "todo")


def mark(st, name, status, bump=None, **extra):
    info = st["stages"].get(name) or {}
    if bump is None:
        bump = 1 if status in ("failed", "partial") else 0       # parked is not a failure
    info.update({"status": status, "at": now(), "attempts": info.get("attempts", 0) + bump})
    if status != "parked":
        for stale in ("reason", "phase", "work"):
            info.pop(stale, None)
    info.update(extra)
    st["stages"][name] = info
    save_state(st)


def invalidate(st, from_stage):
    keep = STAGES[:STAGES.index(from_stage)]
    st["stages"] = {k: v for k, v in st["stages"].items() if k in keep}


def is_parked(st):
    return any(stage_status(st, s) == "parked" for s in STAGES)


def is_quarantined(st, args):
    if args.retry_quarantined:
        return False
    return any(stage_status(st, s) == "failed" and (st["stages"][s].get("attempts", 0) >= MAX_ATTEMPTS)
               for s in STAGES)


def volatile_files(key):
    d = m1.DATA
    return ([os.path.join(d, n % key) for n in
             ("pubchem_%s_summaries.json", "pubchem_%s_assays.csv", "pubchem_%s_summaries.jsonl",
              "pubchem_%s_lookups.json", "chembl_%s_assays_raw.json", "chembl_%s_assays.csv")]
            + glob.glob(os.path.join(LABELS, safe(key) + "_labeled.*.partial.jsonl")))


def drop_volatile_cache(key):
    """Forget the lists that can grow over time. Downloaded assay descriptions (immutable) are kept."""
    for path in volatile_files(key):
        if os.path.exists(path):
            os.remove(path)
    pages = os.path.join(m1.DATA, "chembl_%s_pages" % key)
    if os.path.isdir(pages):
        shutil.rmtree(pages)


# ---------------------------------------------------------------- C1: fetch (chunked, budgeted)

def geneid_of(entry):
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GeneID":
            try:
                return int(ref["id"])
            except (KeyError, ValueError):
                pass
    return None


def pubchem_aids(kind, value):
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/%s/%s/aids/JSON"
           % (kind, urllib.parse.quote(str(value))))
    try:
        data, _ = get_json_h(url)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return []
        raise
    return data.get("IdentifierList", {}).get("AID", [])


def pc_summary_batch(batch, poison):
    """Summaries for a list of AIDs. A batch that breaks is halved until the bad record is isolated and skipped."""
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/summary/JSON" % ",".join(str(a) for a in batch))
    try:
        data, _ = get_json_h(url, tries=1 if len(batch) > 1 else 2, pause=1.0)     # a failing batch is halved, not nursed
        return list(data["AssaySummaries"]["AssaySummary"])
    except (urllib.error.HTTPError, ValueError, KeyError, TypeError) as error:
        if len(batch) == 1:
            poison.append(batch[0])
            return []
        error_kind = error
    except (urllib.error.URLError, OSError) as error:               # timeouts, resets
        if len(batch) == 1:
            tag = "s%d" % batch[0]
            bad_add(tag)
            if bad_count(tag) >= 2:                                 # this single record times out every time
                poison.append(batch[0])
                return []
            raise
        error_kind = error
    mid = len(batch) // 2
    return pc_summary_batch(batch[:mid], poison) + pc_summary_batch(batch[mid:], poison)


def load_summaries(key, aids):
    items, done, poison = {}, set(), set()
    old = read_json(os.path.join(m1.DATA, "pubchem_%s_summaries.json" % key))
    if isinstance(old, dict):
        items = {int(k): v for k, v in old.items()}
        done = set(items)
    for obj in read_lines(os.path.join(m1.DATA, "pubchem_%s_summaries.jsonl" % key)):
        if "_done" in obj:
            done.update(obj["_done"])
        elif "_poison" in obj:
            poison.update(obj["_poison"])
        elif "AID" in obj:
            try:
                items[int(obj["AID"])] = obj
            except (TypeError, ValueError):
                pass
    return items, done, poison


def pubchem_skip_description(summary):
    """ChEMBL copies (read from ChEMBL instead) and siRNA/genome-wide screens are dropped later; do not pay for them."""
    if not summary:
        return False
    if str(summary.get("SourceName", "")).lower() == "chembl":
        return True
    return bool(m2.SIRNA.search(m1.to_text(summary.get("Name")) + " " + m1.to_text(summary.get("Description"))[:600]))


def fetch_pubchem(row, key, entry, args, mode, st):
    acc, gene = row["accession"], row["gene"]
    folder = m1.DATA
    csv_path = os.path.join(folder, "pubchem_%s_assays.csv" % key)
    if args.pubchem_local:                              # PubChem comes from the downloaded zips: no network at all
        rows = read_csv(os.path.join(CACHE, "pubchem_bulk", "csv", key + ".csv"))
        write_csv(csv_path, rows, PUBCHEM_COLUMNS)
        return {"pubchem_assays": len(rows), "notes": []}
    if os.path.exists(csv_path):
        return {"pubchem_assays": len(read_csv(csv_path)), "notes": []}
    lookups_path = os.path.join(folder, "pubchem_%s_lookups.json" % key)
    lookups = read_json(lookups_path, {}) or {}
    kinds = ([("genesymbol", gene)] if gene else []) + [("accession", acc)] + \
            ([("geneid", geneid_of(entry))] if geneid_of(entry) else [])
    notes, last = [], None
    for kind, value in kinds:
        if lookups.get(kind) is not None:
            continue
        try:
            lookups[kind] = sorted(set(pubchem_aids(kind, value)))
        except (urllib.error.URLError, OSError, ValueError) as error:
            lookups[kind] = None
            last = error
            notes.append("lookup by %s failed (%s)" % (kind, type(error).__name__))
        write_json(lookups_path, lookups)
        BUDGET.sleep(args.pause * PACE[0])
    if all(lookups.get(k) is None for k, _ in kinds):
        raise last or OSError("no PubChem lookup worked")
    aids = sorted(set().union(*[set(v) for v in lookups.values() if v]))
    if not aids:
        return {"pubchem_assays": 0, "notes": notes}
    extra = {}
    if mode["cap"] and len(aids) > mode["cap"]:
        if mode["name"] == "light":
            WORK[0] = len(aids) // SUMMARY_BATCH + 40
            raise ParkedCap(len(aids), mode["cap"])
        extra["truncated"] = "PubChem assays sampled down to %d of %d (heavy cap)" % (mode["cap"], len(aids))
        aids = m1.pick_sample(aids, mode["cap"])
    items, done, poison = load_summaries(key, aids)
    missing = [a for a in aids if a not in done]
    for i in range(0, len(missing), SUMMARY_BATCH):
        WORK[0] = (len(missing) - i) // SUMMARY_BATCH + 40
        BUDGET.check()
        batch = missing[i:i + SUMMARY_BATCH]
        bad = []
        got = pc_summary_batch(batch, bad)
        append_lines(os.path.join(folder, "pubchem_%s_summaries.jsonl" % key),
                     got + [{"_done": batch}] + ([{"_poison": bad}] if bad else []))
        for item in got:
            try:
                items[int(item["AID"])] = item
            except (KeyError, TypeError, ValueError):
                pass
        done.update(batch)
        poison.update(bad)
        BUDGET.sleep(args.pause * PACE[0])

    cands = [a for a in aids if not pubchem_skip_description(items.get(a))]
    limit = args.pubchem_desc_limit
    chosen = [] if limit == 0 else m1.pick_sample(cands, limit if limit and limit > 0 else None)
    consecutive = 0
    for n, aid in enumerate(chosen):
        WORK[0] = len(chosen) - n
        path = os.path.join(folder, "pubchem_raw", "AID_%d.json" % aid)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            continue
        tag = "d%d" % aid
        if bad_count(tag) >= 2:
            continue
        BUDGET.check()
        try:
            cached_json("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%d/description/JSON" % aid, path, tries=2)
            consecutive = 0
        except (urllib.error.URLError, OSError, ValueError) as error:
            bad_add(tag)
            consecutive += 1
            if consecutive >= 5 and not isinstance(error, urllib.error.HTTPError):
                raise
            if consecutive >= 12:
                raise
        BUDGET.sleep(args.pause * PACE[0])

    rows = []
    for aid in aids:
        s = items.get(aid, {})
        raw = read_json(os.path.join(folder, "pubchem_raw", "AID_%d.json" % aid))
        f = m1.parse_description(raw) if raw else {}
        title = m1.to_text(s.get("Name")) or f.get("name", "")
        found_by = "|".join(k for k, v in lookups.items() if v and aid in v)
        rows.append({"aid": aid, "source": m1.to_text(s.get("SourceName")), "title": title, "found_by": found_by,
                     "summary_text": m1.clip(s.get("Description"), 1500),
                     "full_description": m1.clip(f.get("description", ""), 1500),
                     "protocol": m1.clip(f.get("protocol", ""), 1500), "targets": f.get("targets", "")})
    write_csv(csv_path, rows, PUBCHEM_COLUMNS)
    if poison:
        extra["skipped_bad_aids"] = sorted(poison)[:20]
    return dict({"pubchem_assays": len(aids), "notes": notes}, **extra)


def lookups_all(lookups):
    return sorted(set().union(*[set(v) for v in lookups.values() if v]))


def chembl_row(a):
    variant = a.get("variant_sequence") or {}
    return {"assay_chembl_id": a.get("assay_chembl_id", ""), "assay_type": a.get("assay_type", ""),
            "description": m1.clip(a.get("description"), 1500), "bao_label": a.get("bao_label", ""),
            "assay_category": a.get("assay_category", ""), "cell_type": a.get("assay_cell_type", ""),
            "organism": a.get("assay_organism", ""),
            "mutation": variant.get("mutation", "") if isinstance(variant, dict) else "",
            "confidence_score": a.get("confidence_score", ""), "document": a.get("document_chembl_id", ""),
            "src_id": a.get("src_id", ""), "src_assay_id": a.get("src_assay_id", "")}


def fetch_chembl(row, key, args, chembl_map, st):
    acc = row["accession"]
    folder = m1.DATA
    csv_path = os.path.join(folder, "chembl_%s_assays.csv" % key)
    if os.path.exists(csv_path):
        return {"chembl_assays": len(read_csv(csv_path)), "notes": []}
    target_id = chembl_map.get(acc)
    if not target_id:                                                # --include-no-chembl: ask directly by accession
        data, _ = get_json_h("https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession=%s"
                             "&target_type=SINGLE%%20PROTEIN" % urllib.parse.quote(acc))
        targets = data.get("targets") or []
        target_id = targets[0]["target_chembl_id"] if targets else None
        if not target_id:
            return {"chembl_assays": 0, "notes": [], "chembl": "no human single-protein target"}
    assays = None
    old = read_json(os.path.join(folder, "chembl_%s_assays_raw.json" % key))       # from version 1
    if isinstance(old, list):
        assays = old
    if assays is None:
        assays, offset, total, skipped = [], 0, None, 0
        pages = os.path.join(folder, "chembl_%s_pages" % key)
        limit = CHEMBL_PAGE
        while True:
            BUDGET.check()
            path = os.path.join(pages, "offset_%07d.json" % offset)
            page = read_json(path)
            if page is None:
                url = ("https://www.ebi.ac.uk/chembl/api/data/assay.json?target_chembl_id=%s&limit=%%d&offset=%d"
                       % (target_id, offset))
                try:
                    page = get_json_h(url % limit, tries=2, pause=1.0)[0]
                except urllib.error.HTTPError as error:
                    if error.code in (500, 502) and limit > 1:      # ChEMBL chokes on one record: smaller pages find it
                        limit = max(1, limit // 4)
                        continue
                    if error.code in (500, 502) and limit == 1:     # that single assay cannot be read: skip it, keep going
                        skipped += 1
                        offset += 1
                        limit = 10
                        continue
                    raise
                write_json(path, page)
                BUDGET.sleep(args.pause * PACE[0])
            meta = page.get("page_meta") or {}
            total = meta.get("total_count", total)
            got = page.get("assays", [])
            assays.extend(got)
            WORK[0] = max(0, ((total or 0) - len(assays)) // CHEMBL_PAGE) + 30
            if not meta.get("next") or not got:
                break
            offset += len(got)
            limit = min(CHEMBL_PAGE, limit * 2)                      # recover to full-size pages after the trouble spot
        shutil.rmtree(pages, ignore_errors=True)
    write_csv(csv_path, [chembl_row(a) for a in assays], CHEMBL_COLUMNS)
    return {"chembl_assays": len(assays), "notes": [], "chembl_skipped_unreadable": skipped} if skipped else \
        {"chembl_assays": len(assays), "notes": []}


def step_fetch(row, key, args, chembl_map, st, mode):
    acc = row["accession"]
    entry = cached_json("https://rest.uniprot.org/uniprotkb/%s.json" % acc,
                        os.path.join(m1.DATA, "uniprot_%s.json" % acc))
    result = {"notes": []}
    if acc in chembl_map or args.include_no_chembl:
        r = fetch_chembl(row, key, args, chembl_map, st)
        result["chembl_assays"] = r.get("chembl_assays", 0)
        if r.get("chembl_skipped_unreadable"):
            result["chembl_skipped_unreadable"] = r["chembl_skipped_unreadable"]
        result["notes"] += r["notes"]
    if not args.skip_pubchem:
        r = fetch_pubchem(row, key, entry, args, mode, st)
        result["notes"] += ["pubchem: " + n for n in r["notes"]]
        for k in ("pubchem_assays", "truncated", "skipped_bad_aids"):
            if k in r:
                result[k] = r[k]
    return result


# ---------------------------------------------------------------- C2: label (sliced, resumable)

def step_label(row, key, args, terms, embedder_box, abstracts, term_hash):
    acc = row["accession"]
    records, dropped = m2.build_records(key, m1.DATA, abstracts)
    if records is None:
        records, dropped = [], []
    if not args.no_original_aids and not args.pubchem_local:      # with local PubChem the originals are already on disk
        missing = sorted({r["pubchem_aid"] for r in records if r["source"] == "ChEMBL"
                          and r["pubchem_aid"].isdigit() and r["text_level"] == "ChEMBL one-line description only"},
                         key=int)[:args.original_aid_limit]
        got = consecutive = 0
        for n, aid in enumerate(missing):
            WORK[0] = len(missing) - n + len(records) // LABEL_SLICE
            path = os.path.join(m1.DATA, "pubchem_raw", "AID_%s.json" % aid)
            if os.path.exists(path) and os.path.getsize(path) > 0:
                continue
            if bad_count("d" + aid) >= 2:
                continue
            BUDGET.check()
            try:
                cached_json("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/description/JSON" % aid, path, tries=2)
                got += 1
                consecutive = 0
            except (urllib.error.URLError, OSError, ValueError) as error:
                bad_add("d" + aid)               # one missing description is not worth failing the protein for
                consecutive += 1
                if consecutive >= 5 and not isinstance(error, urllib.error.HTTPError):
                    raise
            BUDGET.sleep(args.pause * PACE[0])
        if got:
            records, dropped = m2.build_records(key, m1.DATA, abstracts)
    partial = os.path.join(LABELS, "%s_labeled.%s.partial.jsonl" % (safe(key), term_hash))
    saved = {o["id"]: o["f"] for o in read_lines(partial) if "id" in o}
    todo = []
    for rec in records:
        if rec["id"] in saved:
            rec.update(saved[rec["id"]])
        else:
            todo.append(rec)
    if todo:
        embedder = embedder_box.get()
        length = m2.protein_length(m1.DATA, acc)
        implied = m2.implied_function(terms, m1.DATA, acc)
        patterns = m2.construct_patterns(key)
        for i in range(0, len(todo), LABEL_SLICE):
            WORK[0] = (len(todo) - i) // LABEL_SLICE + 1
            BUDGET.check()
            chunk = todo[i:i + LABEL_SLICE]
            with contextlib.redirect_stdout(io.StringIO()):
                m2.label_records(chunk, terms, embedder, patterns, length, implied)
            append_lines(partial, [{"id": r["id"], "f": {k: r.get(k, "") for k in LABEL_FIELDS}} for r in chunk])
    write_csv(os.path.join(LABELS, key + "_labeled.csv"), records, LABEL_COLUMNS)
    for path in glob.glob(os.path.join(LABELS, safe(key) + "_labeled.*.partial.jsonl")):
        os.remove(path)
    return {"records": len(records), "dropped": len(dropped)}


# ---------------------------------------------------------------- C3: summary (with campaign de-duplication)

def merge_effect(effects):
    """Direction of a campaign from the directions of its records. A blank direction counts as unspecified
    (before this fix a blank made the campaign 'mixed', which put it in both the inhibition and activation counts)."""
    e = set(effects) - {"unspecified", ""}
    if e == {"inhibition"}:
        return "inhibition"
    if e == {"activation"}:
        return "activation"
    return "mixed" if e else "unspecified"


def pubchem_originals(records):
    """PubChem AID -> campaign key of the PubChem original, for every PubChem record in the list.
    Used with record_campaign() so a ChEMBL copy of a PubChem deposit is counted with its original.
    (One shared helper for dedupe_campaigns, plate_formats.py, conditions.py and equipment_scan.py.)"""
    original = {}
    for r in records:
        if (r.get("source") or "").startswith("PubChem/") and r.get("pubchem_aid"):
            original[str(r["pubchem_aid"]).strip()] = r.get("campaign") or r["id"]
    return original


def record_campaign(record, originals):
    """The campaign a labeled record belongs to (ChEMBL copies of PubChem deposits join the PubChem original)."""
    camp = record.get("campaign") or record["id"]
    aid = str(record.get("pubchem_aid") or "").strip()
    if record.get("source") == "ChEMBL" and aid in originals:
        return originals[aid]
    return camp


PATENT_KIND = re.compile(r"(?<=\d)[A-Z]\d?$")       # kind code at the end of a patent number: B2, A1, E ...


def canonical_campaign(key, aliases=None):
    """One spelling per campaign, so copies of the same study from different databases count once.
    * aliases (from bindingdb_label.py) map a ChEMBL document to the patent BindingDB says it is;
    * patent numbers lose punctuation and the kind code: PubChem's 'US9023767B2' = BindingDB's 'US9023767'."""
    key = (aliases or {}).get(key, key)
    if key.startswith("patent:"):
        number = re.sub(r"[^A-Za-z0-9]", "", key[len("patent:"):]).upper()
        key = "patent:" + PATENT_KIND.sub("", number)
    return key


def dedupe_campaigns(labeled):
    """A ChEMBL assay that is a copy of a PubChem deposit and the PubChem original are ONE campaign."""
    original = pubchem_originals(labeled)
    camps = collections.OrderedDict()
    merged = 0
    for r in labeled:
        camp = record_campaign(r, original)
        if camp != (r.get("campaign") or r["id"]):
            merged += 1
        c = camps.setdefault(camp, {"f": set(), "t": set(), "e": [], "x": set(), "r": set()})
        c["f"].update(split(r.get("function_labels")))
        c["t"].update(split(r.get("technology_labels")))
        c["e"].append(r.get("effect", ""))
        c["x"].add(r.get("context", ""))
        c["r"].add(r.get("role", ""))
    campaigns = [{"c": k, "f": sorted(v["f"]), "t": sorted(v["t"]), "e": merge_effect(v["e"]),
                  "x": "cell-based" if "cell-based" in v["x"] else "cell-free", "r": sorted(v["r"])}
                 for k, v in camps.items()]
    return campaigns, merged


def scope_of(row, key, counts):
    gene = (row.get("gene") or "").strip()
    return gene.upper() if gene and counts.get(gene, 0) <= 1 else row["accession"]


def step_summary(row, key, classes_rows, scope):
    acc = row["accession"]
    entry = read_json(os.path.join(m1.DATA, "uniprot_%s.json" % acc), {})
    info = m4.summarize_entry(entry) if entry else {"ecs": [], "keywords": [], "go": [], "features": []}
    found = m4.classify(info, classes_rows) if entry else {}
    labeled = read_csv(os.path.join(LABELS, key + "_labeled.csv"))
    campaigns, merged = dedupe_campaigns(labeled)
    summary = {"acc": acc, "gene": key, "scope": scope, "classes": [c for c in found if c != "generic"],
               "n_records": len(labeled), "n_campaigns": len(campaigns), "n_copies_merged": merged,
               "campaigns": campaigns, "labeled_hash": file_hash(os.path.join(LABELS, key + "_labeled.csv")),
               "built": now()}
    write_json(os.path.join(SUMMARY, safe(acc) + ".json"), summary)
    return {"campaigns": len(campaigns), "classes": summary["classes"]}


class EmbedderBox:
    """Loads the meaning model once, only if a protein actually needs labeling."""

    def __init__(self, enabled):
        self.enabled, self.obj, self.tried = enabled, None, False

    def get(self):
        if self.enabled and not self.tried:
            self.tried = True
            try:
                t0 = time.time()
                self.obj = m2.Embedder(LABELS)
                log("  meaning model ready in %.1f s" % (time.time() - t0))
            except Exception as error:
                log("  meaning model not available (%s: %s); using synonyms only" % (type(error).__name__, error))
        return self.obj

    def save(self):
        if self.obj is not None:
            try:
                self.obj.save()
            except Exception:
                pass


# ---------------------------------------------------------------- C: one protein, under a budget

def process_protein(row, key, args, ctx, mode):
    """Returns (status, outcome) with status one of done / parked / failed / quarantined."""
    acc = row["accession"]
    st = load_state(acc)
    term_hash = ctx["term_hash"]
    if st.get("term_hash") and st["term_hash"] != term_hash and stage_done(st, "label"):
        invalidate(st, "label")                         # term table edited: relabel, no downloads
    cutoff = ctx["stale_before"]
    fetched = (st["stages"].get("fetch") or {}).get("at", "")
    if (acc in ctx["refresh"] or (row["gene"] or "").upper() in ctx["refresh"]
            or (cutoff is not None and stage_done(st, "fetch") and fetched <= cutoff)):
        drop_volatile_cache(key)
        st = {"acc": acc, "stages": {}, "v": STATE_VERSION}
    if (stage_done(st, "fetch") and (acc in ctx["chembl_map"] or args.include_no_chembl)
            and not os.path.exists(os.path.join(m1.DATA, "chembl_%s_assays.csv" % key))):
        invalidate(st, "fetch")                         # version 1 missed ChEMBL for this protein: repair (cheap)
    outcome = []
    BUDGET.start_protein(mode["seconds"])
    try:
        for name in STAGES:
            if stage_done(st, name):
                continue
            info = st["stages"].get(name) or {}
            if info.get("status") == "failed" and info.get("attempts", 0) >= MAX_ATTEMPTS and not args.retry_quarantined:
                outcome.append("%s: quarantined after %d failed runs (%s)" % (name, info["attempts"], info.get("error", "")[:60]))
                return "quarantined", outcome
            if info.get("status") == "failed" and args.no_retry:
                outcome.append("%s: skipped (failed earlier)" % name)
                return "failed", outcome
            started = time.time()
            WORK[0] = 0
            try:
                if name == "fetch":
                    if info.get("status") == "partial":         # retry the lists, keep every cached chunk
                        for path in (os.path.join(m1.DATA, "pubchem_%s_assays.csv" % key),
                                     os.path.join(m1.DATA, "chembl_%s_assays.csv" % key)):
                            if os.path.exists(path):
                                os.remove(path)
                    res = step_fetch(row, key, args, ctx["chembl_map"], st, mode)
                    status = "partial" if res["notes"] and info.get("attempts", 0) + 1 < MAX_ATTEMPTS else "done"
                    mark(st, name, status, **res)
                    if status == "partial":
                        outcome.append("fetch partial (%s)" % res["notes"][0])
                        return "failed", outcome
                elif name == "label":
                    res = step_label(row, key, args, ctx["terms"], ctx["embedder"], ctx["abstracts"], term_hash)
                    st["term_hash"] = term_hash
                    mark(st, name, "done", **res)
                else:
                    res = step_summary(row, key, ctx["classes_rows"], ctx["scopes"].get(acc, (row["gene"] or acc).upper()))
                    mark(st, name, "done", **{"campaigns": res["campaigns"]})
                outcome.append("%s ok %.1fs" % (name, time.time() - started))
                ctx["net_failures"] = 0
            except (BudgetExceeded, ParkedCap) as stop:
                reason = stop.reason if isinstance(stop, BudgetExceeded) else str(stop)
                mark(st, name, "parked", reason=reason, phase=mode["name"], work=int(WORK[0]))
                ctx["embedder"].save()
                if isinstance(stop, BudgetExceeded) and stop.kind == "global":
                    outcome.append("%s paused (run time limit)" % name)
                    raise
                outcome.append("%s PARKED (%s)" % (name, reason))
                return "parked", outcome
            except NetworkDown:
                mark(st, name, "failed", bump=0, error="host kept failing")      # an outage is not the protein's fault
                ctx["net_failed"].append((row, key))
                raise
            except (urllib.error.URLError, urllib.error.HTTPError, OSError) as error:
                mark(st, name, "failed", bump=0, error="%s: %s" % (type(error).__name__, str(error)[:150]))
                ctx["net_failed"].append((row, key))
                ctx["net_failures"] += 1
                outcome.append("%s FAILED (%s %s)" % (name, type(error).__name__, str(error)[:50]))
                if ctx["net_failures"] >= args.stop_after_failures:
                    raise NetworkDown("%d failures in a row" % ctx["net_failures"])
                return "failed", outcome
            except Exception as error:                      # a bug or odd record: note it, move on
                mark(st, name, "failed", error="%s: %s %s" % (type(error).__name__, str(error)[:150], getattr(error, "url", "") or ""))
                outcome.append("%s FAILED (%s: %s)" % (name, type(error).__name__, str(error)[:80]))
                return "failed", outcome
        return "done", outcome
    finally:
        BUDGET.end_protein()


# ---------------------------------------------------------------- D: aggregate

def norm(text):
    return re.sub(r"[^a-z0-9]+", " ", str(text).lower()).strip()


def summary_scope(s):
    if s.get("scope"):
        return s["scope"]
    if s["gene"].endswith("_" + s["acc"]):                  # version 1 summaries: duplicate gene names
        return s["acc"]
    return s["gene"].upper()


def merge_campaign_lists(a, b):
    """Union of two campaign lists; a campaign with the same key (e.g. the same patent) is merged, not counted twice."""
    out = collections.OrderedDict()
    for c in list(a) + list(b):
        m = out.get(c["c"])
        if m is None:
            out[c["c"]] = {"c": c["c"], "f": set(c["f"]), "t": set(c["t"]), "e": [c["e"]], "x": {c["x"]}, "r": set(c.get("r", [])),
                           "p": set(c.get("p", [])), "k": dict(c.get("k") or {})}
        else:
            m["f"].update(c["f"]); m["t"].update(c["t"]); m["e"].append(c["e"]); m["x"].add(c["x"]); m["r"].update(c.get("r", []))
            m["p"].update(c.get("p", []))
            for key, val in (c.get("k") or {}).items():
                m["k"].setdefault(key, val)
    return [{"c": k, "f": sorted(v["f"]), "t": sorted(v["t"]), "e": merge_effect(v["e"]),
             "x": "cell-based" if "cell-based" in v["x"] else "cell-free", "r": sorted(v["r"]), "p": sorted(v["p"]), "k": v["k"]}
            for k, v in out.items()]


EFFECTS = ("inhibition", "activation")

CONDS_DIR = os.path.join(CACHE, "conditions")        # written by conditions.py (optional)
CONDITION_COLUMNS = ["stat_id", "entry_id", "scope_type", "scope_id", "parameter", "kind", "n_campaigns_using",
                     "n_campaigns_stating", "median", "p25", "p75", "min", "max", "top_values", "source",
                     "retrieved_on", "method_version", "notes"]
COND_MIN_N = {"gene": 10, "function_class": 5}       # campaigns that must state a parameter before a row is written

PLATES_DIR = os.path.join(CACHE, "plates")          # written by plate_formats.py (optional)
PLATE_RE = re.compile(r"(?<![\d.,])(?:(?P<n>96|384|1536|3456|6144)[\s\-]*(?:well|wells|microtit\w*|micro-?plate\w*|mtp|plates?|format)\b"
                      r"|(?P<w>six|6|12|24|48)[\s\-]+well)", re.I)
PLATE_BEFORE = re.compile(r"(source|compound|mother|stock|master|daughter|dilution|intermediate|library)\b[^.;]{0,25}$", re.I)
PLATE_AFTER = re.compile(r"[\s\-]*(?:well\s*)?(?:\w+\s+){0,2}?(?:source|compound|mother|stock|master|daughter|dilution)\s+plate", re.I)
PLATE_NAME = {"96": "96", "384": "384", "1536": "1536", "3456": "1536", "6144": "1536"}


def plates_in(text):
    """Assay plate formats a text states: '96', '384', '1536' (1536 and denser) and 'low' (6 to 48 wells).
    Mentions of compound / source / stock plates are ignored."""
    found = set()
    text = str(text or "")
    for m in PLATE_RE.finditer(text):
        if PLATE_BEFORE.search(text[max(0, m.start() - 40):m.start()]) or PLATE_AFTER.match(text[m.end():m.end() + 40]):
            continue
        found.add(PLATE_NAME[m.group("n")] if m.group("n") else "low")
    return sorted(found)


# ---------------------------------------------------------------- assay conditions stated in the text (pH, temperature, buffer ...)

_U = r"[µμu]"                                              # micro sign variants
COND_NUMERIC = ["ph", "incubation_min", "dmso_pct", "mgcl2_mm", "nacl_mm", "atp_um", "bsa_pct", "volume_ul"]
COND_CATEGORICAL = ["temperature", "buffer", "detergent", "reducing_agent"]
COND_LABEL = {"ph": "pH", "incubation_min": "incubation time (min)", "dmso_pct": "DMSO (% final)", "mgcl2_mm": "MgCl2 (mM)",
              "nacl_mm": "NaCl (mM)", "atp_um": "ATP (uM)", "bsa_pct": "BSA (%)", "volume_ul": "assay volume (uL)",
              "temperature": "temperature", "buffer": "buffer", "detergent": "detergent", "reducing_agent": "reducing agent"}
_PH = re.compile(r"\bpH\s*(?:of|=|:|~|about|approximately)?\s*(\d{1,2}(?:\.\d{1,2})?)", re.I)
_TEMP_C = re.compile(r"(?<![\d.])(\d{2})\s*(?:°|º|˚|degrees?)\s*(?:celsius|c)\b|(?<![\d.])(\d{2})\s*c\b(?!\w)", re.I)
_RT = re.compile(r"room temp|ambient temp|\brt\b|\bambient\b", re.I)
_TIME = re.compile(r"(?:incubat\w*|reaction|react\w*|for|after)\W{0,40}?(\d+(?:\.\d+)?)\s*(min(?:ute)?s?|h(?:ou)?rs?|h)\b", re.I)
_DMSO1 = re.compile(r"(\d+(?:\.\d+)?)\s*%\s*(?:\(?v\s*/\s*v\)?\s*)?(?:final\s+)?DMSO", re.I)
_DMSO2 = re.compile(r"DMSO\W{0,25}?(\d+(?:\.\d+)?)\s*%", re.I)
_MG = re.compile(r"(\d+(?:\.\d+)?)\s*mM\s*(?:of\s*)?MgCl\s*2|MgCl\s*2\W{0,8}?(\d+(?:\.\d+)?)\s*mM", re.I)
_NA = re.compile(r"(\d+(?:\.\d+)?)\s*mM\s*(?:of\s*)?NaCl|NaCl\W{0,8}?(\d+(?:\.\d+)?)\s*mM", re.I)
_ATP = re.compile(r"(\d+(?:\.\d+)?)\s*(" + _U + r"M|mM|nM)\s*(?:of\s*)?ATP\b|\bATP\W{0,12}?(\d+(?:\.\d+)?)\s*(" + _U + r"M|mM|nM)", re.I)
_BSA1 = re.compile(r"(\d*\.?\d+)\s*%\s*(?:\(?w\s*/\s*v\)?\s*)?(?:BSA|bovine serum albumin)", re.I)
_BSA2 = re.compile(r"(\d*\.?\d+)\s*mg\s*/\s*m[lL]\s*(?:BSA|bovine serum albumin)", re.I)
_VOL = re.compile(r"(?:final|total|assay|reaction)\s+(?:assay\s+|reaction\s+)?volume\W{0,25}?(\d+(?:\.\d+)?)\s*" + _U + r"[lL]\b", re.I)
_BUFFERS = [("HEPES", r"hepes"), ("Tris", r"\btris\b|tris-?hcl"), ("phosphate / PBS", r"\bpbs\b|phosphate buffer|sodium phosphate|potassium phosphate"),
            ("MOPS", r"\bmops\b"), ("MES", r"\bmes\b"), ("bicine", r"bicine"), ("HBSS", r"hbss|hank'?s"), ("citrate", r"citrate buffer|sodium citrate"),
            ("glycine", r"glycine"), ("CHES", r"\bches\b"), ("acetate", r"acetate buffer|sodium acetate")]
_DETERGENTS = [("Tween-20", r"tween[- ]?20|polysorbate"), ("Triton X-100", r"triton"), ("Brij-35", r"brij"), ("CHAPS", r"chaps"),
               ("NP-40 / Nonidet", r"nonidet|np-?40"), ("Pluronic", r"pluronic"), ("CHAPSO / octyl glucoside", r"chapso|octyl[- ]glucoside"), ("SDS", r"\bsds\b")]
_REDUCERS = [("DTT", r"\bdtt\b|dithiothreitol"), ("TCEP", r"\btcep\b"), ("2-mercaptoethanol", r"mercaptoethanol|\b2-?me\b|\bbme\b"), ("glutathione", r"glutathione|\bgsh\b")]
_BUF_C = [(n, re.compile(rx, re.I)) for n, rx in _BUFFERS]
_DET_C = [(n, re.compile(rx, re.I)) for n, rx in _DETERGENTS]
_RED_C = [(n, re.compile(rx, re.I)) for n, rx in _REDUCERS]


def _first(rx, text, group=1):
    m = rx.search(text)
    return float(m.group(group) or m.group(group + 1)) if m else None


def conditions_in(text):
    """Assay conditions a text states (first stated value; categories as lists).
    Returns {parameter: value}; numbers are floats. Only plausible values are kept."""
    text = str(text or "")
    out = {}
    for m in _PH.finditer(text):
        v = float(m.group(1))
        if 3.0 <= v <= 11.0:
            out["ph"] = v
            break
    for m in _TEMP_C.finditer(text):
        v = int(m.group(1) or m.group(2))
        if 4 <= v <= 60:
            out["temperature"] = ["%d C" % v]
            break
    else:
        if _RT.search(text):
            out["temperature"] = ["room temperature"]
    for m in _TIME.finditer(text):
        v = float(m.group(1)) * (60 if m.group(2).lower().startswith("h") else 1)
        if 1 <= v <= 1440:
            out["incubation_min"] = v
            break
    for rx in (_DMSO1, _DMSO2):
        vals = [float(m.group(1)) for m in rx.finditer(text)]
        vals = [v for v in vals if 0 < v <= 10]            # 100 % DMSO is the compound stock, not the assay
        if vals:
            out["dmso_pct"] = vals[0]
            break
    v = _first(_MG, text)
    if v is not None and 0 < v <= 100:
        out["mgcl2_mm"] = v
    v = _first(_NA, text)
    if v is not None and 0 < v <= 1000:
        out["nacl_mm"] = v
    m = _ATP.search(text)
    if m:
        num = float(m.group(1) or m.group(3))
        unit = (m.group(2) or m.group(4)).lower().replace("µ", "u").replace("μ", "u")
        um = num * 1000 if unit == "mm" else num / 1000 if unit == "nm" else num
        if 0 < um <= 10000:
            out["atp_um"] = um
    m = _BSA1.search(text)
    if m:
        out["bsa_pct"] = float(m.group(1))
    else:
        m = _BSA2.search(text)
        if m:
            out["bsa_pct"] = float(m.group(1)) / 10.0
    if "bsa_pct" in out and not 0 < out["bsa_pct"] <= 10:
        del out["bsa_pct"]
    v = _first(_VOL, text)
    if v is not None and 0 < v <= 1000:
        out["volume_ul"] = v
    for key, table in (("buffer", _BUF_C), ("detergent", _DET_C), ("reducing_agent", _RED_C)):
        found = [n for n, rx in table if rx.search(text)]
        if found:
            out[key] = found
    return out


def quantile(sorted_values, q):
    if not sorted_values:
        return None
    pos = (len(sorted_values) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(sorted_values) - 1)
    return sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * (pos - lo)


def effect_subset(campaigns, effect):
    """Campaigns that speak to the wanted direction.
    inhibition: everything not labeled activation (binding and unlabeled assays are normally run to find inhibitors/binders);
    activation: only campaigns labeled activation, or mixed (both directions)."""
    if effect == "activation":
        return [c for c in campaigns if c["e"] in ("activation", "mixed")]
    return [c for c in campaigns if c["e"] != "activation"]


def aggregate(args, assays, classes_rows, stamp):
    summaries = [read_json(p) for p in
                 (os.path.join(SUMMARY, f) for f in sorted(os.listdir(SUMMARY)) if f.endswith(".json"))] \
        if os.path.isdir(SUMMARY) else []
    summaries = [s for s in summaries if s]
    n_plates = 0
    for s in summaries:                                         # plate formats found by plate_formats.py
        pl = read_json(os.path.join(PLATES_DIR, safe(s["acc"]) + ".json"))
        if pl:
            n_plates += 1
            for c in s["campaigns"]:
                c["p"] = sorted(set(c.get("p", [])) | set(pl["campaigns"].get(c["c"], [])))
    n_conds = 0
    for s in summaries:                                         # assay conditions found by conditions.py
        cd = read_json(os.path.join(CONDS_DIR, safe(s["acc"]) + ".json"))
        if cd:
            n_conds += 1
            for c in s["campaigns"]:
                c["k"] = cd["campaigns"].get(c["c"]) or c.get("k") or {}
    # One key per study, then add BindingDB. This runs AFTER plates and conditions are attached, because those
    # files use the original keys. Copies of one study (PubChem patent deposit, BindingDB patent row, ChEMBL
    # document of the same patent) become one campaign; see canonical_campaign().
    use_bdb = not getattr(args, "no_bindingdb", False) and os.path.isdir(SUMMARY_BDB)
    n_bdb = 0
    for s in summaries:
        extra = (read_json(os.path.join(SUMMARY_BDB, safe(s["acc"]) + ".json")) if use_bdb else None) or {}
        aliases = extra.get("aliases") or {}
        own = [dict(c, c=canonical_campaign(c["c"], aliases)) for c in s["campaigns"]]
        other = [dict(c, c=canonical_campaign(c["c"])) for c in extra.get("campaigns") or []]
        n_bdb += bool(other)
        s["campaigns"] = merge_campaign_lists(own, other)
        s["n_campaigns"] = len(s["campaigns"])
    entries =[(a, set(split(a.get("evidence_labels") or a.get("evidence_function")))) for a in assays]
    entries = [(a, labels) for a, labels in entries if labels]
    version = "bd-" + stamp.replace("-", "")
    source = "build_database.py (ChEMBL + PubChem%s, labeled by m2)" % (" + BindingDB" if n_bdb else "")
    today = datetime.now().strftime("%Y-%m-%d")
    stat_rows, effect_rows, cond_rows = [], [], []

    def rows_for(scope_type, scope_id, campaigns_all, note_extra="", classes=()):
        out = []
        for effect in EFFECTS:
            out += rows_for_effect(scope_type, scope_id, effect_subset(campaigns_all, effect), effect, note_extra, classes)
        return out

    def using_for(scope_type, scope_id, campaigns, classes):
        classifiable = [c for c in campaigns if c["f"]]
        for a, wanted in entries:
            entry_classes = set(split(a.get("classes")))
            if scope_type == "function_class" and scope_id not in entry_classes:
                continue
            if scope_type == "gene" and not (entry_classes & (set(classes) | {"generic"})):
                continue                          # a kinase assay is not "used" for a protein that is not a kinase
            scope_text = (a.get("scope") or "").lower()
            want_cell = "cell-based" in scope_text and "free" not in scope_text
            # a cell-based assay entry is only counted from cell-based campaigns, and a cell-free one from cell-free
            using = [c for c in classifiable if wanted & set(c["f"]) and (c.get("x") == "cell-based") == want_cell]
            if using:
                yield a, wanted, using

    def cond_rows_for(scope_type, scope_id, campaigns, classes):
        out = []
        for a, wanted, using in using_for(scope_type, scope_id, campaigns, classes):
            ks = [c.get("k") or {} for c in using]
            common = [source, today, version, "auto-generated from stated values in assay texts; see notes"]
            for param in COND_NUMERIC:
                vals = sorted(k[param] for k in ks if param in k)
                if len(vals) < COND_MIN_N[scope_type]:
                    continue
                out.append(dict(zip(CONDITION_COLUMNS, ["%s-%s-%s" % (scope_id, a["entry_id"], param), a["entry_id"], scope_type,
                    scope_id, param, "numeric", len(using), len(vals), round(quantile(vals, .5), 3), round(quantile(vals, .25), 3),
                    round(quantile(vals, .75), 3), vals[0], vals[-1], ""] + common)))
            for param in COND_CATEGORICAL:
                stating = [k[param] for k in ks if k.get(param)]
                if len(stating) < COND_MIN_N[scope_type]:
                    continue
                cnt = collections.Counter(v for lst in stating for v in set(lst))
                top = "; ".join("%s|%.3f" % (v, n / len(stating)) for v, n in cnt.most_common(5))
                out.append(dict(zip(CONDITION_COLUMNS, ["%s-%s-%s" % (scope_id, a["entry_id"], param), a["entry_id"], scope_type,
                    scope_id, param, "categorical", len(using), len(stating), "", "", "", "", "", top] + common)))
        return out

    def rows_for_effect(scope_type, scope_id, campaigns, effect, note_extra, classes):
        out = []
        total = len(campaigns)
        classifiable = [c for c in campaigns if c["f"]]
        for a, wanted, using in using_for(scope_type, scope_id, campaigns, classes):
            base = [a["entry_id"], scope_type, scope_id, "any", effect]
            tail = [source, "functions: " + "; ".join(sorted(wanted)), today, version, AUTO_ROW_NOTE + note_extra]
            sid = "%s-%s-%s" % (scope_id, a["entry_id"], effect[:3])
            out.append(dict(zip(USAGE_COLUMNS, [sid] + base + [""] + [total, len(classifiable), len(using), "", ""] + tail)))
            tech = collections.Counter(t for c in using for t in set(c["t"]))
            n_tech = sum(1 for c in using if c["t"])
            for t, n in tech.most_common():
                if n < args.min_readout and len(tech) >= 5:
                    continue                      # one-off technology labels are mostly text noise
                out.append(dict(zip(USAGE_COLUMNS, [sid + "-" + re.sub(r"\W+", "", t)[:12]] + base + [t]
                                    + [total, n_tech, n, "", ""] + tail)))
            plate_c = [c for c in using if c.get("p")]           # campaigns whose text states a plate format
            for fmt, n in collections.Counter(p for c in plate_c for p in set(c["p"])).most_common():
                out.append(dict(zip(USAGE_COLUMNS, [sid + "-plate" + fmt] + base + ["plate:" + fmt]
                                    + [total, len(plate_c), n, "", ""] + tail)))
        return out

    class_campaigns = collections.defaultdict(list)    # class -> campaigns of all its proteins (merged below)
    gene_count = collections.defaultdict(set)
    for s in summaries:
        scope = summary_scope(s)
        stat_rows += rows_for("gene", scope, s["campaigns"], classes=s["classes"])
        cond_rows += cond_rows_for("gene", scope, s["campaigns"], s["classes"])
        counts = collections.Counter(c["e"] for c in s["campaigns"])
        effect_rows.append({"gene": scope, "accession": s["acc"], "classes": "; ".join(s["classes"]),
                            "n_campaigns": s["n_campaigns"], "inhibition": counts["inhibition"],
                            "activation": counts["activation"], "mixed": counts["mixed"],
                            "unspecified": counts["unspecified"],
                            "no_function_label": sum(1 for c in s["campaigns"] if not c["f"])})
        for cls in s["classes"]:
            gene_count[cls].add(scope)
            class_campaigns[cls].extend(s["campaigns"])
    class_merged = {}
    for cls, camps in class_campaigns.items():
        merged = merge_campaign_lists(camps, [])        # a campaign shared by several proteins counts once
        class_merged[cls] = merged
        stat_rows += rows_for("function_class", cls, merged, "; counted across %d proteins" % len(gene_count[cls]))
        cond_rows += cond_rows_for("function_class", cls, merged, ())

    seen, unique = {}, []                           # de-duplicate on the key the database uses
    for r in stat_rows:
        key = (r["entry_id"], r["scope_type"], r["scope_id"], r["stage"], r["effect"], r["readout"])
        if key not in seen:
            seen[key] = True
            unique.append(r)
    write_csv(os.path.join(OUTDIR, "usage_stats.generated.csv"), unique, USAGE_COLUMNS)
    write_csv(os.path.join(OUTDIR, "condition_stats.generated.csv"), cond_rows, CONDITION_COLUMNS)
    write_csv(os.path.join(OUTDIR, "effect_by_gene.csv"), effect_rows,
              ["gene", "accession", "classes", "n_campaigns", "inhibition", "activation", "mixed", "unspecified",
               "no_function_label"])
    class_rows, gap_rows = [], []
    for cls, merged in sorted(class_merged.items()):
        counts = collections.Counter(c["e"] for c in merged)
        class_rows.append({"class_id": cls, "n_proteins": len(gene_count[cls]), "n_campaigns": len(merged),
                           "inhibition": counts["inhibition"], "activation": counts["activation"],
                           "mixed": counts["mixed"], "unspecified": counts["unspecified"]})
        gap_rows.append({"class_id": cls, "n_proteins": len(gene_count[cls]), "n_campaigns": len(merged),
                         "no_function_label": sum(1 for c in merged if not c["f"]),
                         "no_technology_label": sum(1 for c in merged if not c["t"])})
    write_csv(os.path.join(OUTDIR, "effect_by_class.csv"), class_rows,
              ["class_id", "n_proteins", "n_campaigns", "inhibition", "activation", "mixed", "unspecified"])
    write_csv(os.path.join(OUTDIR, "function_gaps.csv"), gap_rows,
              ["class_id", "n_proteins", "n_campaigns", "no_function_label", "no_technology_label"])

    # readout suggestions: technologies seen often for a class that no playbook entry for that class lists
    readouts_by_class = collections.defaultdict(list)
    for a in assays:
        for c in split(a.get("classes")):
            readouts_by_class[c].append((a["entry_id"], norm(a.get("readouts", ""))))
    suggestions = []
    for cls, camps in class_merged.items():
        pair = collections.Counter()
        for c in camps:
            for f in c["f"]:
                for t in c["t"]:
                    pair[(f, t)] += 1
        for (f, t), n in pair.most_common():
            if n < args.suggest_min:
                continue
            parts = [norm(p) for p in re.split(r"[/(),]", t) if len(norm(p)) > 2] + [norm(t)]
            # generic entries (thermal shift, SPR, ITC ...) apply to every class, so they count as covering it too
            listed = readouts_by_class.get(cls, []) + readouts_by_class.get("generic", [])
            covered = sorted({e for e, text in listed if any(p in text for p in parts)})
            suggestions.append({"class_id": cls, "function_label": f, "technology": t, "n_campaigns": n,
                                "n_proteins_in_class": len(gene_count[cls]),
                                "playbook_entries_listing_it": "; ".join(covered),
                                "suggestion": "covered" if covered else "NOT in any playbook readout list for this class"})
    write_csv(os.path.join(OUTDIR, "readout_suggestions.csv"), suggestions,
              ["class_id", "function_label", "technology", "n_campaigns", "n_proteins_in_class",
               "playbook_entries_listing_it", "suggestion"])
    no_labels = [a["entry_id"] for a in assays if not (a.get("evidence_labels") or a.get("evidence_function"))]
    return {"proteins": len(summaries), "with_bindingdb": n_bdb, "with_plates": n_plates, "with_conditions": n_conds, "condition_rows": len(cond_rows), "rows": len(unique),
            "suggestions": sum(1 for s in suggestions if s["suggestion"] != "covered"),
            "entries_without_evidence_labels": no_labels}


# ---------------------------------------------------------------- progress report (parked / quarantined)

def progress_report(names):
    rows, counts = [], collections.Counter()
    if not os.path.isdir(STATE):
        return rows, counts
    for f in sorted(os.listdir(STATE)):
        if not f.endswith(".json"):
            continue
        st = read_json(os.path.join(STATE, f), {}) or {}
        acc = st.get("acc", f[:-5])
        statuses = [stage_status(st, s) for s in STAGES]
        for s, v in zip(STAGES, statuses):
            counts["%s:%s" % (s, v)] += 1
        if all(v == "done" for v in statuses):
            continue
        bad = next((s for s, v in zip(STAGES, statuses) if v != "done"), "")
        info = st["stages"].get(bad) or {}
        status = info.get("status", "todo")
        if status == "failed" and info.get("attempts", 0) >= MAX_ATTEMPTS:
            status = "quarantined"
        rows.append({"accession": acc, "gene": names.get(acc, ""), "stage": bad, "status": status,
                     "reason": info.get("reason") or info.get("error", ""), "work_left_estimate": info.get("work", ""),
                     "attempts": info.get("attempts", 0), "last_update": info.get("at", "")})
    return rows, counts


def write_progress(names):
    rows, counts = progress_report(names)
    write_csv(os.path.join(OUTDIR, "parked_proteins.csv"), rows,
              ["accession", "gene", "stage", "status", "reason", "work_left_estimate", "attempts", "last_update"])
    return rows, counts


# ---------------------------------------------------------------- E: merge

def summarized_scopes():
    """(scope_type, scope_id) of every protein and function class the aggregate step looked at.
    An automatic row for one of these scopes that the latest build no longer produces is out of date."""
    out = set()
    if os.path.isdir(SUMMARY):
        for f in os.listdir(SUMMARY):
            s = read_json(os.path.join(SUMMARY, f)) if f.endswith(".json") else None
            if s:
                out.add(("gene", summary_scope(s)))
                out.update(("function_class", c) for c in s.get("classes", []))
    return out


def merge_into_data(replace_m2):
    gen = read_csv(os.path.join(OUTDIR, "usage_stats.generated.csv"))
    path = os.path.join(DATA, "usage_stats.csv")
    existing = read_csv(path)

    def key_of(r):
        return (r["entry_id"], r["scope_type"], r["scope_id"], r["stage"], r.get("effect") or "any", r["readout"])

    def is_auto(r):
        v = r.get("method_version", "")
        return v.startswith(AUTO_PREFIXES) or (replace_m2 and v == "m2-labels")

    index = {key_of(r): i for i, r in enumerate(existing)}
    gen_keys = {key_of(g) for g in gen}
    # scopes this build covered: those with new rows, plus summarized proteins that now produce no rows at all
    # (before, their old automatic rows stayed forever; 83 such rows were found in Nov 2026, e.g. LGALS8)
    regenerated = {(g["scope_type"], g["scope_id"]) for g in gen} | summarized_scopes()
    added = replaced = kept = 0
    conflicts = []
    for g in gen:
        k = key_of(g)
        if k not in index:
            existing.append(g)
            index[k] = len(existing) - 1
            added += 1
        elif is_auto(existing[index[k]]):
            existing[index[k]] = g
            replaced += 1
        else:
            kept += 1
            conflicts.append("%s / %s %s" % (g["entry_id"], g["scope_id"], g["readout"] or "(whole assay)"))
    # earlier automatic rows that this build no longer produces (inside scopes it did regenerate) are removed
    pruned = 0
    final = []
    for r in existing:
        if is_auto(r) and key_of(r) not in gen_keys and (r["scope_type"], r["scope_id"]) in regenerated:
            pruned += 1
            continue
        final.append(r)
    for r in final:
        r["effect"] = r.get("effect") or "any"          # rows written before the effect column existed
    backup = os.path.join(DATA, "backup")
    os.makedirs(backup, exist_ok=True)
    if os.path.exists(path):
        shutil.copy(path, os.path.join(backup, "usage_stats_%s.csv" % datetime.now().strftime("%Y%m%d_%H%M%S")))
    write_csv(path, final, USAGE_COLUMNS)
    cgen = os.path.join(OUTDIR, "condition_stats.generated.csv")
    cpath = os.path.join(DATA, "condition_stats.csv")
    if os.path.exists(cgen) and read_csv(cgen):               # the whole table is generated; keep a backup of the last one
        if os.path.exists(cpath):
            shutil.copy(cpath, os.path.join(backup, "condition_stats_%s.csv" % datetime.now().strftime("%Y%m%d_%H%M%S")))
        shutil.copy(cgen, cpath)
    return added, replaced, kept, conflicts, pruned


# ---------------------------------------------------------------- main

def build_parser():
    ap = argparse.ArgumentParser(description="One-hit database build (resumable, updateable, time-boxed per protein).")
    ap.add_argument("--genes", help="comma-separated gene symbols or accessions (default: all)")
    ap.add_argument("--genes-file", help="text file with one gene symbol or accession per line")
    ap.add_argument("--limit", type=int, help="work on at most this many proteins in the light pass")
    ap.add_argument("--max-minutes", type=float, help="hard limit for the whole run (also inside a protein); rerun to continue")
    ap.add_argument("--light-minutes", type=float, default=6.0, help="time budget per protein in pass 1 (default 6)")
    ap.add_argument("--light-cap", type=int, default=2500, help="PubChem assays above which a protein is parked in pass 1")
    ap.add_argument("--heavy-minutes", type=float, default=40.0, help="time budget per parked protein in pass 2 (default 40)")
    ap.add_argument("--heavy-cap", type=int, default=15000, help="PubChem assays kept per protein in pass 2 (0 = no cap)")
    ap.add_argument("--light-only", action="store_true", help="skip pass 2 (parked proteins wait for a later run)")
    ap.add_argument("--heavy-only", action="store_true", help="skip pass 1; work on parked proteins only")
    ap.add_argument("--retry-quarantined", action="store_true", help="try proteins that failed 3 runs in a row again")
    ap.add_argument("--include-no-chembl", action="store_true", help="also try proteins with no ChEMBL target")
    ap.add_argument("--pubchem-local", action="store_true",
                    help="use the downloaded PubChem zips (run pubchem_plan.py and pubchem_local.py first); no PubChem calls")
    ap.add_argument("--skip-pubchem", action="store_true", help="ChEMBL only (much faster)")
    ap.add_argument("--pubchem-desc-limit", type=int, default=40,
                    help="full PubChem descriptions to download per protein (default 40; 0 = titles only)")
    ap.add_argument("--no-original-aids", action="store_true", help="skip downloading PubChem originals of ChEMBL copies")
    ap.add_argument("--original-aid-limit", type=int, default=60)
    ap.add_argument("--rules-only", action="store_true", help="labels from term-table synonyms only (no meaning model)")
    ap.add_argument("--refresh", help="comma-separated genes to download again from scratch")
    ap.add_argument("--refresh-older-than", type=int, metavar="DAYS", help="download again anything fetched more than DAYS ago")
    ap.add_argument("--refresh-list", action="store_true", help="download the UniProt and ChEMBL lists again")
    ap.add_argument("--no-retry", action="store_true", help="do not retry steps that failed earlier")
    ap.add_argument("--stop-after-failures", type=int, default=8,
                    help="this many network failures in a row = rest and try again (default 8)")
    ap.add_argument("--timeout", type=float, default=45, help="seconds to wait for one request (default 45; try 90 on a slow link)")
    ap.add_argument("--cooldown-seconds", type=float, default=120,
                    help="first rest after a run of failures; doubles each time, up to 15 min (default 120)")
    ap.add_argument("--max-cooldowns", type=int, default=6,
                    help="rests allowed per run before stopping (default 6)")
    ap.add_argument("--pause", type=float, default=0.4, help="seconds between calls (be polite to the services)")
    ap.add_argument("--min-readout", type=int, default=2,
                    help="drop a technology row used by fewer campaigns than this (only when 5 or more technologies are seen)")
    ap.add_argument("--suggest-min", type=int, default=5, help="minimum campaigns for a readout suggestion")
    ap.add_argument("--aggregate-only", action="store_true", help="skip downloading/labeling; just rebuild the output tables")
    ap.add_argument("--no-aggregate", action="store_true")
    ap.add_argument("--merge", action="store_true", help="add new rows to data/usage_stats.csv (backup is made)")
    ap.add_argument("--replace-m2", action="store_true", help="with --merge, also replace earlier m2-labels rows")
    ap.add_argument("--no-bindingdb", action="store_true", help="ignore campaigns from bindingdb_label.py")
    ap.add_argument("--status", action="store_true", help="only show progress so far")
    return ap


def needs_work(st, p, ctx):
    if not all(stage_done(st, s) for s in STAGES):
        return True
    if st.get("term_hash") != ctx["term_hash"]:
        return True
    if p["accession"] in ctx["refresh"] or (p["gene"] or "").upper() in ctx["refresh"]:
        return True
    return ctx["stale_before"] is not None and (st["stages"]["fetch"].get("at", "") <= ctx["stale_before"])


def run_pass(label, items, args, ctx, mode, tally, limit=None):
    """Work through items. Returns 'time' if the global time limit stopped it, else ''."""
    for n, (p, key) in enumerate(items, 1):
        if BUDGET.global_deadline and time.time() >= BUDGET.global_deadline:
            log("  run time limit reached; rerun the same command to continue")
            return "time"
        if limit and tally["worked"] >= limit:
            return ""
        tally["worked"] += 1
        t0 = time.time()
        while True:
            try:
                status, outcome = process_protein(p, key, args, ctx, mode)
                break
            except BudgetExceeded:
                log("  run time limit reached; rerun the same command to continue")
                return "time"
            except NetworkDown as down:
                if tally["cooldowns"] >= args.max_cooldowns:
                    raise
                tally["cooldowns"] += 1
                wait = min(900, args.cooldown_seconds * 2 ** (tally["cooldowns"] - 1))
                log("  connection trouble (%s). Resting %d s, then trying again (rest %d of %d)"
                    % (down, wait, tally["cooldowns"], args.max_cooldowns))
                ctx["net_failures"] = 0
                BREAKER.fails.clear()
                BREAKER.opens.clear()
                BREAKER.until.clear()
                PACE[0] = max(PACE[0], 3.0)
                try:
                    BUDGET.sleep(wait)
                except BudgetExceeded:
                    log("  run time limit reached; rerun the same command to continue")
                    return "time"
        tally[status] += 1
        log("  [%s %d/%d] %-10s %-8s %s  (%.1fs)" % (label, n, len(items), p["gene"] or "-", p["accession"],
                                                    "; ".join(outcome), time.time() - t0))
        if tally["worked"] % 25 == 0:
            ctx["embedder"].save()
    return ""


def main(argv=None):
    args = build_parser().parse_args(argv)
    for folder in (CACHE, STATE, SUMMARY, LABELS, LOGS, OUTDIR, os.path.join(CACHE, "m1_data")):
        os.makedirs(folder, exist_ok=True)
    stamp = datetime.now().strftime("%Y-%m-%d")
    log.handle = open(os.path.join(LOGS, "run_%s.log" % datetime.now().strftime("%Y%m%d_%H%M%S")), "w", encoding="utf-8")
    route_m1_through_http_get()
    m1.DATA = os.path.join(CACHE, "m1_data")
    PAUSE[0] = args.pause
    TIMEOUT[0] = args.timeout
    m2.say = lambda *a, **k: None
    try:
        signal.signal(signal.SIGTERM, lambda *a: (_ for _ in ()).throw(KeyboardInterrupt()))
    except (ValueError, OSError, AttributeError):
        pass
    if args.pubchem_local and not os.path.isdir(os.path.join(CACHE, "pubchem_bulk", "csv")):
        sys.exit("--pubchem-local needs build_cache/pubchem_bulk/csv: run pubchem_plan.py and pubchem_local.py first.")
    assays = read_csv(os.path.join(DATA, "assays.csv"))
    classes_rows = m4.read_table(os.path.join(DATA, "function_classes.csv"))
    log("Build started %s (build_database.py v2)" % now())
    names = {p["accession"]: p["gene"] for p in read_csv(os.path.join(CACHE, "uniprot_human_reviewed.csv"))}

    if args.status:
        rows, counts = progress_report(names)
        log("Proteins with a state file: %d" % len(os.listdir(STATE)))
        for k, v in sorted(counts.items()):
            log("  %-18s %d" % (k, v))
        by = collections.Counter(r["status"] for r in rows)
        log("Not finished: %s" % (", ".join("%d %s" % (v, k) for k, v in by.items()) or "none"))
        for r in sorted(rows, key=lambda r: (r["status"], str(r["work_left_estimate"]).zfill(8)))[:15]:
            log("  %-10s %-8s %-11s %s" % (r["gene"] or "-", r["accession"], r["status"], r["reason"][:70]))
        write_progress(names)
        return

    summary_line = ""
    if not args.aggregate_only:
        BUDGET.global_deadline = time.time() + args.max_minutes * 60 if args.max_minutes else None
        log("Step A: protein list")
        proteins = load_protein_list(args.refresh_list)
        names = {p["accession"]: p["gene"] for p in proteins}
        log("  %d reviewed human proteins" % len(proteins))
        log("Step B: ChEMBL targets")
        chembl_map = load_chembl_targets(args.refresh_list)
        log("  %d accessions with a ChEMBL single-protein target" % len(chembl_map))

        wanted = None
        gnames = [g.strip().upper() for g in (args.genes or "").split(",") if g.strip()]
        if args.genes_file:
            with open(args.genes_file, encoding="utf-8") as handle:
                gnames += [l.strip().upper() for l in handle if l.strip() and not l.startswith("#")]
        if gnames:
            wanted = set(gnames)
        counts = collections.Counter(p["gene"] for p in proteins if p["gene"])
        todo, scopes = [], {}
        for p in proteins:
            if wanted is not None and p["accession"] not in wanted and p["gene"].upper() not in wanted:
                continue
            if wanted is None and not args.include_no_chembl and p["accession"] not in chembl_map:
                continue
            key = safe(p["gene"]) if p["gene"] and counts[p["gene"]] == 1 else safe((p["gene"] or "NA") + "_" + p["accession"])
            todo.append((p, key))
            scopes[p["accession"]] = scope_of(p, key, counts)
        if wanted:
            missing = wanted - {p["accession"] for p, _ in todo} - {p["gene"].upper() for p, _ in todo}
            if missing:
                log("  not found in the human reviewed list: %s" % ", ".join(sorted(missing)))
        log("  %d proteins selected" % len(todo))

        added = ensure_terms()
        if added:
            log("  term table: added %d new terms for the wider playbook (see notes column)" % added)
        terms = m2.load_terms()
        embedder = EmbedderBox(not args.rules_only)
        ctx = {"terms": terms, "term_hash": file_hash(m2.TERM_TABLE) + ("r" if args.rules_only else "m"),
               "embedder": embedder, "abstracts": m2.load_abstracts(m1.DATA), "chembl_map": chembl_map,
               "classes_rows": classes_rows, "net_failures": 0, "net_failed": [], "scopes": scopes,
               "refresh": {g.strip().upper() for g in (args.refresh or "").split(",") if g.strip()},
               "stale_before": (datetime.now() - timedelta(days=args.refresh_older_than)).replace(microsecond=0).isoformat()
               if args.refresh_older_than is not None else None}
        tally = collections.Counter()
        skipped = 0
        stopped = ""
        try:
            if not args.heavy_only:
                light, quarantined = [], 0
                for p, key in todo:
                    st = load_state(p["accession"])
                    if not needs_work(st, p, ctx):
                        skipped += 1
                    elif is_parked(st) and not (p["accession"] in ctx["refresh"] or p["gene"].upper() in ctx["refresh"]):
                        continue                                  # waits for the heavy pass
                    elif is_quarantined(st, args):
                        quarantined += 1
                    else:
                        light.append((p, key))
                log("Pass 1 (light): %d proteins to do, %d already finished, %d quarantined (--retry-quarantined to try again)"
                    % (len(light), skipped, quarantined))
                log("  each protein gets %g min and up to %d PubChem assays; bigger ones are parked, not failed"
                    % (args.light_minutes, args.light_cap))
                mode = {"name": "light", "seconds": args.light_minutes * 60, "cap": args.light_cap}
                stopped = run_pass("light", light, args, ctx, mode, tally, args.limit)
            if not args.light_only and not stopped:
                heavy = []
                for p, key in todo:
                    st = load_state(p["accession"])
                    if needs_work(st, p, ctx) and is_parked(st) and not is_quarantined(st, args):
                        stage = next(s for s in STAGES if stage_status(st, s) == "parked")
                        heavy.append(((st["stages"][stage].get("work") or 0), p, key))
                heavy.sort(key=lambda x: x[0])                    # smallest remaining job first: most proteins finish
                if heavy:
                    log("Pass 2 (heavy): %d parked proteins, smallest remaining job first; %g min and %s assays each"
                        % (len(heavy), args.heavy_minutes, args.heavy_cap or "unlimited"))
                    mode = {"name": "heavy", "seconds": args.heavy_minutes * 60, "cap": args.heavy_cap}
                    stopped = run_pass("heavy", [(p, k) for _, p, k in heavy], args, ctx, mode, tally)
            if not stopped and ctx["net_failed"] and not args.light_only:
                again, seen = [], set()
                for row, key in ctx["net_failed"]:
                    if row["accession"] not in seen and not stage_done(load_state(row["accession"]), "summary"):
                        seen.add(row["accession"])
                        again.append((row, key))
                if again:
                    log("Retry pass: %d proteins that hit network trouble earlier in this run" % len(again))
                    ctx["net_failed"] = []
                    stopped = run_pass("retry", again, args, ctx,
                                       {"name": "light", "seconds": args.light_minutes * 60, "cap": args.light_cap}, tally)
        except NetworkDown as error:
            log("\nSTOPPED: the connection kept failing even after resting (%s). Nothing is lost: rerun the same command later." % error)
            summary_line = "stopped early (network)"
        except KeyboardInterrupt:
            log("\nStopped by you. Nothing is lost: rerun the same command to continue.")
            summary_line = "stopped by user"
        finally:
            embedder.save()
        if stopped == "time":
            summary_line = "stopped at the time limit"
        log("Step C summary: %d proteins finished this run, %d already finished before, %d parked (big or slow; "
            "continue next run), %d failed (retried next run), %d quarantined"
            % (tally["done"], skipped, tally["parked"], tally["failed"], tally["quarantined"]))
        rows, _ = write_progress(names)
        waiting = [r for r in rows if r["status"] == "parked"]
        if waiting:
            log("  %d proteins still parked (list: build_output/parked_proteins.csv). Run again; "
                "add --heavy-minutes 90 to give them more time each." % len(waiting))
        BUDGET.global_deadline = None

    if not args.no_aggregate:
        log("Step D: aggregate (usage rows, effect counts, readout suggestions)")
        res = aggregate(args, assays, classes_rows, stamp)
        log("  %d proteins summarized (%d with BindingDB campaigns) -> %d de-duplicated usage rows (build_output/usage_stats.generated.csv)"
            % (res["proteins"], res["with_bindingdb"], res["rows"]))
        log("  plate formats read for %d proteins%s" % (res["with_plates"], "" if res["with_plates"] else
            " (run scripts/plate_formats.py first to add them)"))
        log("  assay conditions read for %d proteins -> %d condition rows (build_output/condition_stats.generated.csv)%s"
            % (res["with_conditions"], res["condition_rows"], "" if res["with_conditions"] else
               " (run scripts/conditions.py first to add them)"))
        log("  effect counts: build_output/effect_by_gene.csv and effect_by_class.csv; gaps: function_gaps.csv")
        log("  readout suggestions not yet in the playbook: %d (build_output/readout_suggestions.csv)" % res["suggestions"])
        if res["entries_without_evidence_labels"]:
            log("  playbook entries with no evidence_labels (cannot get usage rows yet): %d" % len(res["entries_without_evidence_labels"]))
    if args.merge:
        log("Step E: merge into data/usage_stats.csv")
        added, replaced, kept, conflicts, pruned = merge_into_data(args.replace_m2)
        log("  %d rows added, %d earlier automatic rows replaced, %d out-of-date automatic rows removed, "
            "%d kept because they are yours" % (added, replaced, pruned, kept))
        for c in conflicts[:10]:
            log("    kept yours: " + c)
        if len(conflicts) > 10:
            log("    ... and %d more" % (len(conflicts) - 10))
    log("Done. Log: %s %s" % (log.handle.name, summary_line))


if __name__ == "__main__":
    main()
