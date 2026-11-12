"""
M2 labeling for the Assay Strategy Planner.

Reads the files M1 saved (m1_data folder) and, for each target:
  1. Merges the PubChem and ChEMBL assays into one list, dropping PubChem's
     copies of ChEMBL assays and its siRNA/RNAi screens (kept in a "dropped" file).
  2. Labels every assay on two axes using the editable term table:
        FUNCTION   = what is measured (ATPase, helicase unwinding, ...)
        TECHNOLOGY = how it is measured (ADP-Glo, SPR, gel, ...)
     First by your synonyms (exact phrases). Then, only for gaps, by meaning
     (a small language model on your computer). Weak matches go to a review file.
  3. Reads protein construct residue ranges out of the assay text, for you to
     spot-check.
  4. Groups assays into campaigns and prints a preview ranking.

Files it creates next to this script:
  term_table.csv            the term table (created once; YOU edit it, it is never overwritten)
  m2_data/                  results and caches
  m2_report.txt             the summary to paste back into the chat

How to run (with the Python inside .venv):
    .\\.venv\\Scripts\\python.exe m2_label.py wrn
    .\\.venv\\Scripts\\python.exe m2_label.py            (WRN and EGFR)
    .\\.venv\\Scripts\\python.exe m2_label.py wrn --rules-only   (skip the language model)

Nothing is downloaded from the databases here (that was M1). The language model
is cached, and every embedding it computes is saved in m2_data, so reruns are fast.
"""

import argparse
import collections
import csv
import hashlib
import glob
import json
import os
import re
import sys
import time
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
TERM_TABLE = os.path.join(HERE, "term_table.csv")

GENES = {"WRN": "Q14191", "EGFR": "P00533"}
MODEL = "sentence-transformers/all-MiniLM-L6-v2"

SEM_MIN = 0.50      # a meaning-based match must score at least this ...
SEM_MARGIN = 0.08   # ... and beat the runner-up by at least this
MAX_CHUNKS = 10     # pieces of text embedded per assay
CHUNK_CHARS = 300
MIN_CONSTRUCT_LEN = 30
REVIEW_CAP = 300
# a paper that merely mentions a crystal structure or mass spec did not necessarily use it as THIS assay's readout
ABSTRACT_SKIP = {"Crystallography / cryo-EM", "Mass spectrometry"}
BINDING_WORDS = re.compile(r"binding|affinity|\bkd\b|\bki\b|dissociation|thermal|melting|adduct|covalent|"
                           r"selectivity|\bcell|cellular|proliferat", re.I)

LOG = []


def say(text=""):
    print(text)
    LOG.append(text)


# ------------------------------------------------------------- term table

TERM_COLUMNS = ["axis", "term", "synonyms", "examples", "suppresses", "applies_to", "notes"]

# Terms that only make sense for assays run in cells. A cell-free assay is never
# given these labels, and a synonym hit on one marks the assay as cell-based.
CELL_ONLY = {
    "Target engagement (cellular)", "Cell proliferation / viability", "Cellular signaling",
    "Apoptosis / DNA damage (cellular)", "siRNA / RNAi screen", "CETSA", "Western blot",
    "Luciferase reporter",
}

# axis, term, synonyms (separated by |), examples (sentences for the meaning match), suppresses
SEED = [
    # ---- FUNCTION: what is measured
    ("function", "ATPase activity",
     "atpase|atp hydrolysis|atp turnover|adp production|adp formation",
     "ATPase activity measured as ADP produced from ATP hydrolysis", ""),
    ("function", "Helicase unwinding",
     "helicase|unwinding|unwind|duplex dna unwinding",
     "helicase unwinding of a forked duplex DNA substrate", ""),
    ("function", "Exonuclease / nuclease activity",
     "exonuclease|nuclease|dna degradation|dna cleavage",
     "exonuclease degradation of a DNA substrate", ""),
    ("function", "Kinase activity",
     "kinase|phosphorylation|phosphorylate|autophosphorylation|phosphotransferase|poly glu tyr|poly glu ala tyr|poly gat|atp incorporation|phosphate incorporation|tyrosine kinase",
     "kinase activity measured by phosphorylation of a peptide substrate", ""),
    ("function", "Protein-DNA binding",
     "dna binding|binding to dna|protein dna binding|electrophoretic mobility shift|emsa|gel shift|fluorescence anisotropy",
     "Binding of the protein to a DNA substrate", ""),
    ("function", "Compound-DNA binding (interference)",
     "thiazole orange|intercalator|intercalation|dna intercalat|compound dna binding|dna binding counterscreen",
     "Counterscreen for compounds that bind the DNA substrate", "Protein-DNA binding"),
    ("function", "Ligand binding affinity",
     "binding affinity|binding constant|dissociation constant|kd|kon|koff|association constant|binding kinetics|competitive binding",
     "Binding affinity of the compound to the purified protein", ""),
    ("function", "Covalent adduct formation",
     "adduct|covalent|intact protein mass|cysteine modification",
     "Covalent adduct formation on the protein measured by mass spectrometry", ""),
    ("function", "Time-dependent / covalent kinetics",
     "kinact|k inact|inactivation rate|saturation rate constant|time dependent inhibition|residence time|jump dilution",
     "Rate of covalent inactivation, kinact over Ki", ""),
    ("function", "Kinome selectivity panel",
     "selectivity interaction|kinomescan|kinome scan|kinase panel|selectivity profiling|kinome profiling",
     "Selectivity profiling across a panel of kinases", ""),
    ("function", "Target engagement (cellular)",
     "target engagement|cetsa|cellular thermal shift|thermal stability in cells",
     "Cellular target engagement measured by thermal stabilization", ""),
    ("function", "Cell proliferation / viability",
     "proliferation|viability|cell growth|cytotoxic|cytotoxicity|celltiter|mtt|cell survival",
     "cell proliferation or viability readout in cultured cells", ""),
    ("function", "Cellular signaling",
     "p erk|perk|p akt|phospho erk|downstream signaling|pathway|erk phosphorylation|erk1 2 phosphorylation|"
     "phosphorylation of erk|phosphorylation of erk1 2|perk1 2|p erk1 2|phospho erk1 2",
     "downstream signaling readout in cultured cells", ""),
    ("function", "Apoptosis / DNA damage (cellular)",
     "apoptosis|caspase|cleaved parp|double strand break|double stands break|p chk2|gamma h2ax",
     "Induction of apoptosis or DNA damage in cells", ""),
    ("function", "siRNA / RNAi screen",
     "sirna|rnai|shrna|crispr|genome wide screen",
     "Genome-wide siRNA screen", ""),
    # ---- TECHNOLOGY: how it is measured
    ("technology", "ADP-Glo luminescence", "adpglo|promega adp",
     "ADP-Glo luminescence readout", ""),
    ("technology", "Kinase-Glo luminescence", "kinaseglo", "Kinase-Glo luminescence readout", ""),
    ("technology", "Transcreener", "transcreener", "Transcreener ADP detection", ""),
    ("technology", "TR-FRET / HTRF",
     "htrf|tr fret|trfret|time resolved fret|time resolved fluorescence|ulight|lance ultra",
     "Time-resolved FRET homogeneous assay", "Fluorogenic / FRET substrate"),
    ("technology", "LanthaScreen", "lanthascreen", "LanthaScreen TR-FRET assay", "TR-FRET / HTRF|Fluorogenic / FRET substrate"),
    ("technology", "AlphaScreen / AlphaLISA", "alphascreen|alphalisa|alpha screen",
     "AlphaScreen proximity assay", ""),
    ("technology", "Fluorescence polarization",
     "fluorescence polarization|fluorescence polarisation|fp assay|anisotropy",
     "Fluorescence polarization binding assay", ""),
    ("technology", "Fluorogenic / FRET substrate",
     "fluor quencher|fluorophore quencher|quencher|quenched|tamra|bhq|cy5|fret assay|fluorogenic|fluorometric|fluorescence based|fluorescence intensity|dequench|z lyte|zlyte|fret based",
     "Fluorophore and quencher labeled DNA substrate, fluorescence increases on unwinding", ""),
    ("technology", "Gel-based (PAGE / EMSA)",
     "gel electrophoresis|polyacrylamide|page analysis|gel based|gel shift|emsa|sds page",
     "Products separated by polyacrylamide gel electrophoresis", ""),
    ("technology", "Radiometric", "radiometric|radioactiv|32p|33p|gamma32p|gamma33p|scintillation|flashplate|filter binding",
     "Radiometric assay with 32P or 33P labeled substrate", ""),
    ("technology", "Mobility shift (Caliper / LabChip)", "mobility shift|caliper|labchip",
     "Microfluidic mobility shift assay", ""),
    ("technology", "KINOMEscan competition binding", "kinomescan|kinome scan|ambit",
     "KINOMEscan active-site competition binding assay", ""),
    ("technology", "ELISA / immunosorbent", "elisa|immunosorbent|immunoassay|tmb",
     "ELISA or immunosorbent based assay", ""),
    ("technology", "SPR / BLI",
     "spr|surface plasmon|biolayer interferometry|bli|biacore",
     "Surface plasmon resonance binding kinetics", ""),
    ("technology", "Thermal shift (DSF)",
     "thermal shift|differential scanning fluorimetry|dsf|sypro|melting temperature|thermofluor",
     "Thermal shift assay with SYPRO Orange dye, change in melting temperature", ""),
    ("technology", "CETSA", "cetsa|cellular thermal shift",
     "Cellular thermal shift assay", ""),
    ("technology", "ITC", "isothermal titration|itc", "Isothermal titration calorimetry", ""),
    ("technology", "Mass spectrometry", "ms analysis|esi ms|lc ms|maldi|intact mass|tof",
     "Mass spectrometry analysis", ""),
    ("technology", "Western blot", "western blot|immunoblot", "Western blotting analysis", ""),
    ("technology", "Luciferase reporter", "luciferase reporter|reporter gene|reporter assay",
     "Luciferase reporter gene assay", ""),
    ("technology", "Crystallography / cryo-EM",
     "x ray crystallography|crystal structure|cryo em|cocrystal",
     "X-ray crystal structure of the complex", ""),
]


def write_seed_table():
    with open(TERM_TABLE, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(TERM_COLUMNS)
        for axis, term, synonyms, examples, suppresses in SEED:
            writer.writerow([axis, term, synonyms, examples, suppresses,
                             "cell-based" if term in CELL_ONLY else "any", ""])


# --------------------------------------------------------- text handling

KIT = [
    (re.compile(r"adp[\s-]?glo(?:\s*[™®]|\(tm\))?(?:\s+kinase)?(?:\s+(?:max|assay))*", re.I),
     " adpglo "),
    (re.compile(r"kinase[\s-]?glo(?:\s*[™®])?(?:\s+(?:plus|max))?(?:\s+assay)?", re.I),
     " kinaseglo "),
]


def norm(text):
    """Lower-case, with kit names fused into one word and separators as spaces."""
    t = str(text or "")
    for regex, replacement in KIT:
        t = regex.sub(replacement, t)
    t = t.lower()
    t = re.sub(r"[^a-z0-9]+", " ", t)   # punctuation and symbols become spaces
    return t.strip()


def clean(text):
    return " ".join(str(text or "").split())


def clip(text, limit):
    text = clean(text)
    return text if len(text) <= limit else text[: limit - 3] + "..."


def chunk(text):
    """Split text into pieces of roughly CHUNK_CHARS for the language model."""
    pieces = [p for p in re.split(r"(?<=[.;])\s+", clean(text)) if p]
    out, current = [], ""
    for piece in pieces:
        if current and len(current) + len(piece) > CHUNK_CHARS:
            out.append(current)
            current = piece
        else:
            current = (current + " " + piece).strip()
    if current:
        out.append(current)
    return out[:MAX_CHUNKS] or [""]


class Term:
    def __init__(self, axis, name, synonyms, examples, suppresses, applies_to="any"):
        self.axis, self.name = axis, name
        self.applies_to = "cell-based" if clean(applies_to).lower() == "cell-based" else "any"
        phrases = [name] + [s for s in synonyms.split("|")]
        seen, self.phrases = set(), []
        for phrase in phrases:
            n = norm(phrase)
            if n and n not in seen:
                seen.add(n)
                self.phrases.append(n)
        alternatives = "|".join(re.escape(p) for p in sorted(self.phrases, key=len, reverse=True))
        self.regex = re.compile(r"(?<![a-z0-9])(?:%s)(?![a-z0-9])" % alternatives)
        # phrases shown to the language model: name, synonyms and example sentences
        self.sem_phrases = [p for p in self.phrases]
        for example in examples.split("|"):
            if clean(example):
                self.sem_phrases.append(clean(example))
        self.suppresses = [clean(s) for s in suppresses.split("|") if clean(s)]


def load_terms():
    if not os.path.exists(TERM_TABLE):
        write_seed_table()
        say("Created the starter term table: %s" % TERM_TABLE)
    terms = []
    with open(TERM_TABLE, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            axis = clean(row.get("axis")).lower()
            name = clean(row.get("term"))
            if axis not in ("function", "technology") or not name:
                continue
            applies = row.get("applies_to") or ("cell-based" if name in CELL_ONLY else "any")
            terms.append(Term(axis, name, row.get("synonyms") or "",
                              row.get("examples") or "", row.get("suppresses") or "", applies))
    return terms


def rule_match(terms, axis, text_norm):
    hits = collections.OrderedDict()
    for term in terms:
        if term.axis != axis:
            continue
        match = term.regex.search(text_norm)
        if match:
            hits[term.name] = match.group(0)
    for term in terms:
        if term.axis == axis and term.name in hits:
            for suppressed in term.suppresses:
                hits.pop(suppressed, None)
    return hits


# ------------------------------------------------------ context and roles

EXPRESSION = re.compile(
    r"(?:expressed|extracted|produced|purified|isolated|obtained|derived)\s+(?:in|from)\s+[^;,.()]{0,70}?\bcells?\b"
    r"|(?:baculovirus|insect)[\s-]*(?:infected\s+)?[\w-]*\s*cells"
    r"|\bsf9\b|\bbl21\b|e(?:scherichia)?\.?\s?coli[^;,.()]{0,30}cells", re.I)
CELLISH = re.compile(
    r"\bcells?\b|cell[\s-]based|cell line|cellular|cetsa|hct[\s-]?116|hek[\s-]?293|\ba431\b|\bhela\b|"
    r"western blot|in vivo|xenograft|proliferation|viability|apoptosis|\bsirna\b|\brnai\b|\bshrna\b", re.I)
SIRNA = re.compile(r"\bsirna\b|\brnai\b|\bshrna\b|crispr|genome[\s-]?wide|whole[\s-]genome", re.I)
COUNTER = re.compile(r"counter[\s-]?screen|selectivity|\bblm\b|off[\s-]target|thiazole orange|interference", re.I)
CONFIRM = re.compile(r"confirmatory|cherry[\s-]?pick|confirmation", re.I)
INHIB = re.compile(r"inhibit|\bic50\b|antagonist|suppress", re.I)
ACTIV = re.compile(r"activat|activator|agonist|\bec50\b|stimulat|enhance", re.I)
MUTANT = re.compile(r"\b[A-Z]\d{2,4}[A-Z]\b|\bmutant|\bmutation|\bdel\d|\d+del\b|\bdel\b|\bins\b|exon ?19", re.I)
PATENT = re.compile(r"\b(?:US\d{7,11}[A-Z]\d?|WO\d{9,12})\b")


ANTIPROLIF = re.compile(r"antiproliferat|cytotoxic|cytostatic|in vivo|xenograft", re.I)
HOSTS = re.compile(r"sf9|sf21|high[\s-]?five|\bhi5\b|bl21|e\.?\s?coli|escherichia|insect|baculovirus", re.I)

# Bump this when a labeling RULE in this file changes (not the term table, which has its own hash).
# build_database.py and bindingdb_label.py put it in the label cache key, so cached labels are redone (no downloads).
RULES_VERSION = "rules3"     # 2026-11-11: membrane preparations and GTPgammaS binding count as cell-free;
                             # 2026-11-12: a protein's EC number can imply its activity term (EC_IMPLIES)

# Assays run on membrane preparations (or [35S]GTPgammaS / Eu-GTP binding, which needs membranes) are cell-free,
# even when the receptor was "expressed in CHO cells" or ChEMBL lists a cell type. Plural "membranes" or an explicit
# preparation word is required, so "membrane potential" or "membrane permeability" (cell assays) do not match.
GTP_BINDING = re.compile(r"\[35s?\]\s*-?\s*gtp|\b35s[\s-]*gtp|gtp\s*-?\s*(?:γ|gamma)\s*-?\s*s?\b|gtpγs|\beu[\s-]?gtp", re.I)
MEMBRANE = re.compile(r"\bmembranes\b|\bmembrane[\s-]+(?:prep\w*|homogenates?|fractions?|suspensions?)\b|"
                      r"\bcell[\s-]+membrane\b", re.I)
MEMBRANE_BAO = {"cell membrane format"}     # not "subcellular format": ChEMBL also uses it for lysates and cell assays
# ... unless the text says the cells were intact or living, or the readout is one only living cells give
WHOLE_CELL = re.compile(r"\b(?:whole|intact|live|living)(?:[\s-]+[\w-]+){0,2}?[\s-]+cells?\b|whole[\s-]cell\b|"
                        r"proliferat|viability|apoptosis|western blot|in vivo|xenograft|\bsirna\b|cetsa", re.I)
# ChEMBL often writes "expressed in HEK293 cell membrane" for ordinary cell assays; these readouts give them away
LIVE_READOUT = re.compile(r"calcium|\bca2\+|\bcamp\b|current|membrane potential|surface|density|reporter|corrector|"
                          r"potentiat|activation of|induced activation|internali[sz]|arrestin|impedance|\bflux\b|"
                          r"uptake|efflux|secretion|fluorescen\w* imaging|integrity|permeab|\bldh\b|cytotox|flipr|"
                          r"patch[\s-]?clamp|electrophysiolog", re.I)


def context_of(text, bao="", cell_type=""):
    """cell-based or cell-free. Protein made in insect or bacterial cells is still cell-free,
    and so is an assay on membranes prepared from cells: [35S]GTPgammaS / Eu-GTP binding, or a text (or ChEMBL
    format) naming membranes with no readout that needs living cells."""
    if not WHOLE_CELL.search(text):
        if GTP_BINDING.search(text):
            return "cell-free"
        if (MEMBRANE.search(text) or bao in MEMBRANE_BAO) and not LIVE_READOUT.search(text):
            return "cell-free"
    stripped = EXPRESSION.sub(" ", text)
    if CELLISH.search(stripped) or ANTIPROLIF.search(stripped):
        return "cell-based"
    if cell_type and not HOSTS.search(cell_type):
        return "cell-based"
    if bao == "cell-based format" and not EXPRESSION.search(text) and not HOSTS.search(text):
        return "cell-based"
    return "cell-free"


# ------------------------------------------------------ construct ranges

DASH = r"(?:to|-|–|—)"
TO_END = re.compile(r"\b[A-Z]?(\d{1,4})\s*to\s*end\b", re.I)   # "696 to end residues"


def construct_patterns(gene):
    return [
        re.compile(r"\b[A-Z]?(\d{1,4})\s*%s\s*[A-Z]?(\d{1,4})\s*(?:residues?|aa|amino acids)\b" % DASH, re.I),
        re.compile(r"\b(?:residues?|aa|amino acids?)\s*[A-Z]?(\d{1,4})\s*%s\s*[A-Z]?(\d{1,4})\b" % DASH, re.I),
        re.compile(r"\b%s\b\s*\(\s*(?:aa\s*)?[A-Z]?(\d{1,4})\s*%s\s*[A-Z]?(\d{1,4})\s*\)" % (re.escape(gene), DASH), re.I),
    ]


def find_construct(text, patterns, length):
    rejected = None
    for pattern in patterns:
        for match in pattern.finditer(text):
            start, end = int(match.group(1)), int(match.group(2))
            if start >= end or end - start + 1 < MIN_CONSTRUCT_LEN:
                rejected = rejected or "range %d-%d too short or reversed" % (start, end)
                continue
            if length and end > length + 5:
                rejected = rejected or "range %d-%d is longer than the protein (%d)" % (start, end, length)
                continue
            lo, hi = max(0, match.start() - 70), min(len(text), match.end() + 70)
            note = "range stated"
            if re.search(r"full[\s-]?length", text, re.I) and length and (end - start + 1) < 0.9 * length:
                note = "range stated (text also says full-length)"
            return {"start": start, "end": end, "note": note, "evidence": clean(text[lo:hi])}
    to_end = TO_END.search(text)
    if to_end and length and 0 < int(to_end.group(1)) < length - MIN_CONSTRUCT_LEN:
        lo, hi = max(0, to_end.start() - 70), min(len(text), to_end.end() + 70)
        return {"start": int(to_end.group(1)), "end": length, "note": "range stated (to the C-terminal end)",
                "evidence": clean(text[lo:hi])}
    if re.search(r"full[\s-]?length", text, re.I):
        return {"start": "", "end": "", "note": "full-length stated, no range", "evidence": ""}
    if re.search(r"helicase domain|kinase domain|catalytic domain|truncat|fragment|\bcore\b", text, re.I):
        return {"start": "", "end": "", "note": "partial (domain named), no range", "evidence": ""}
    return {"start": "", "end": "", "note": rejected or "not stated", "evidence": ""}


# ------------------------------------------------------------ read M1 data


def read_csv(path):
    if not os.path.exists(path):
        return None
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def protein_length(m1_dir, accession):
    path = os.path.join(m1_dir, "uniprot_%s.json" % accession)
    try:
        with open(path, encoding="utf-8") as handle:
            return int(json.load(handle)["sequence"]["length"])
    except (OSError, KeyError, ValueError, TypeError):
        return None


def raw_pubchem_text(path):
    """Name, description and protocol from a saved PubChem description file."""
    try:
        with open(path, encoding="utf-8") as handle:
            descr = json.load(handle)["PC_AssayContainer"][0]["assay"]["descr"]
    except (OSError, KeyError, IndexError, TypeError, ValueError):
        return ""
    parts = [descr.get("name"), descr.get("description"), descr.get("protocol"), descr.get("comment")]
    return clean(" ".join(clean(p) for p in parts if p))


def build_records(gene, m1_dir, abstracts=None):
    """Merge ChEMBL and PubChem rows; return (kept, dropped)."""
    chembl = read_csv(os.path.join(m1_dir, "chembl_%s_assays.csv" % gene))
    pubchem = read_csv(os.path.join(m1_dir, "pubchem_%s_assays.csv" % gene))
    if chembl is None and pubchem is None:
        return None, None
    kept, dropped = [], []
    seen_ids = set()
    for row in chembl or []:
        rid = "chembl:" + row["assay_chembl_id"]
        if rid in seen_ids:
            continue
        seen_ids.add(rid)
        campaign = "doc:" + row["document"] if row.get("document") else rid
        one_liner = clean(row.get("description"))
        text, level = one_liner, "ChEMBL one-line description only"
        aid = clean(row.get("src_assay_id")) if str(row.get("src_id")) == "7" else ""
        if aid.isdigit():
            extra = raw_pubchem_text(os.path.join(m1_dir, "pubchem_raw", "AID_%s.json" % aid))
            if extra:
                text = one_liner + " || PubChem AID %s: %s" % (aid, clip(extra, 3000))
                level = "ChEMBL + original PubChem text"
        kept.append({
            "id": rid, "source": "ChEMBL", "campaign": campaign, "text_level": level,
            "text": text, "title": one_liner,
            "chembl_type": row.get("assay_type", ""), "bao": row.get("bao_label", ""),
            "cell_type": row.get("cell_type", ""), "mutation": row.get("mutation", ""),
            "pubchem_aid": row.get("src_assay_id", "") if str(row.get("src_id")) == "7" else "",
            "abstract": (abstracts or {}).get(row.get("document", ""), ""),
        })
    for row in pubchem or []:
        rid = "pubchem:" + row["aid"]
        source = clean(row.get("source"))
        title = clean(row.get("title"))
        text = clean(" ".join([title, row.get("full_description", ""), row.get("protocol", "")]))
        level = "PubChem full text" if clean(row.get("full_description")) else "PubChem title only"
        record = {"id": rid, "source": "PubChem/" + source, "text": text, "title": title,
                  "text_level": level,
                  "chembl_type": "", "bao": "", "mutation": "", "pubchem_aid": row["aid"]}
        if source.lower() == "chembl":
            record["reason"] = "copy of a ChEMBL assay"
            dropped.append(record)
            continue
        if SIRNA.search(title + " " + clean(row.get("summary_text"))[:600]):
            record["reason"] = "siRNA / RNAi / genome-wide screen"
            dropped.append(record)
            continue
        patent = PATENT.search(title)
        record["campaign"] = ("patent:" + patent.group(0)) if patent else "pc:" + row["aid"]
        kept.append(record)
    return kept, dropped



# ------------------------------------------------------ M2b: abstracts and implied function

def load_abstracts(m1_dir):
    """Paper title + abstract per ChEMBL document id, saved by m2b_fetch_abstracts.py."""
    found = {}
    for path in glob.glob(os.path.join(m1_dir, "chembl_docs", "*.json")):
        try:
            with open(path, encoding="utf-8") as handle:
                doc = json.load(handle)
        except (OSError, ValueError):
            continue
        text = clean(" ".join(clean(doc.get(k)) for k in ("title", "abstract") if doc.get(k)))
        if text:
            found[doc.get("document_chembl_id") or os.path.basename(path)[:-5]] = text
    return found


EC_IMPLIES = [            # (EC prefix, function term). The first matching prefix decides for each EC number.
    ("3.6.5.", None),            # GTPases: "inhibition of KRAS" says nothing about which activity was read
    ("7.", None),                # translocases (ABC transporters, ion pumps): measured as transport, not turnover
    ("3.4.19.12", "Deubiquitinase activity"),
    ("3.4.", "Protease activity"),
    ("2.7.10.", "Kinase activity"), ("2.7.11.", "Kinase activity"), ("2.7.12.", "Kinase activity"),
    # lipid kinases (PI3K, PI4K, PIP kinases, DGK, SPHK): same assay platforms as protein kinases
    ("2.7.1.137", "Kinase activity"), ("2.7.1.153", "Kinase activity"), ("2.7.1.154", "Kinase activity"),
    ("2.7.1.67", "Kinase activity"), ("2.7.1.68", "Kinase activity"), ("2.7.1.149", "Kinase activity"),
    ("2.7.1.150", "Kinase activity"), ("2.7.1.107", "Kinase activity"), ("2.7.1.91", "Kinase activity"),
    ("2.7.1.138", "Kinase activity"), ("2.7.1.158", "Kinase activity"),
    # histone demethylases (KDM1, JmjC KDM2-7) and histone methyltransferases (SET, DOT1L, PRMT)
    ("1.14.11.27", "Histone demethylase activity"), ("1.14.11.65", "Histone demethylase activity"),
    ("1.14.11.66", "Histone demethylase activity"), ("1.14.11.67", "Histone demethylase activity"),
    ("1.14.11.68", "Histone demethylase activity"), ("1.14.11.69", "Histone demethylase activity"),
    ("1.14.99.66", "Histone demethylase activity"),
    ("2.1.1.43", "Histone methyltransferase activity"), ("2.1.1.319", "Histone methyltransferase activity"),
    ("2.1.1.320", "Histone methyltransferase activity"), ("2.1.1.321", "Histone methyltransferase activity"),
    ("2.1.1.322", "Histone methyltransferase activity"), ("2.1.1.354", "Histone methyltransferase activity"),
    ("2.1.1.355", "Histone methyltransferase activity"), ("2.1.1.356", "Histone methyltransferase activity"),
    ("2.1.1.357", "Histone methyltransferase activity"), ("2.1.1.359", "Histone methyltransferase activity"),
    ("2.1.1.360", "Histone methyltransferase activity"), ("2.1.1.361", "Histone methyltransferase activity"),
    ("2.1.1.362", "Histone methyltransferase activity"), ("2.1.1.364", "Histone methyltransferase activity"),
    ("2.1.1.367", "Histone methyltransferase activity"), ("2.1.1.369", "Histone methyltransferase activity"),
    ("2.1.1.370", "Histone methyltransferase activity"), ("2.1.1.371", "Histone methyltransferase activity"),
    ("3.1.3.16", "Phosphatase activity"), ("3.1.3.48", "Phosphatase activity"),
    ("3.1.4.17", "Phosphodiesterase activity"), ("3.1.4.35", "Phosphodiesterase activity"),
    ("3.1.4.53", "Phosphodiesterase activity"),
    ("1.14.14.", "Cytochrome P450 activity"), ("1.14.15.", "Cytochrome P450 activity"),
    ("2.3.1.48", "Histone acetyltransferase activity"),
    ("2.1.1.37", "DNA / RNA methyltransferase activity"), ("2.1.1.348", "DNA / RNA methyltransferase activity"),
    ("3.5.1.98", "Histone deacetylase activity"),
    ("2.4.2.30", "ADP-ribosyltransferase activity"), ("2.4.2.31", "ADP-ribosyltransferase activity"),
    ("2.3.2.", "Ubiquitination activity"),
    ("5.6.2.", "Topoisomerase activity"),
    ("3.1.11.", "Exonuclease / nuclease activity"), ("3.1.13.", "Exonuclease / nuclease activity"),
    ("3.1.14.", "Exonuclease / nuclease activity"), ("3.1.15.", "Exonuclease / nuclease activity"),
    ("3.1.16.", "Exonuclease / nuclease activity"), ("3.1.21.", "Exonuclease / nuclease activity"),
    ("3.1.22.", "Exonuclease / nuclease activity"), ("3.1.25.", "Exonuclease / nuclease activity"),
    ("3.1.26.", "Exonuclease / nuclease activity"), ("3.1.27.", "Exonuclease / nuclease activity"),
    ("3.1.30.", "Exonuclease / nuclease activity"), ("3.1.31.", "Exonuclease / nuclease activity"),
    ("", "Enzyme activity"),     # any other EC number
]


def ec_numbers(entry):
    """EC numbers in a UniProt JSON entry (recommended and alternative names, and domains/components)."""
    out = []

    def walk(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k == "ecNumbers" and isinstance(v, list):
                    out.extend(x.get("value", "") for x in v if isinstance(x, dict))
                else:
                    walk(v)
        elif isinstance(obj, list):
            for item in obj:
                walk(item)
    walk(entry.get("proteinDescription", {}))
    return sorted(set(e for e in out if e))


def ec_implied(ecs, term_names):
    """The one function term the EC numbers point to, or None (no EC, an excluded EC, or two different terms)."""
    found = set()
    for ec in ecs:
        for prefix, term in EC_IMPLIES:
            exact = prefix and prefix[-1].isdigit()          # "2.1.1.37" must not match "2.1.1.370"
            if (ec == prefix) if exact else ec.startswith(prefix):
                found.add(term)
                break
    if None in found:
        return None
    found = {t for t in found if t in term_names}
    return found.pop() if len(found) == 1 else None


def implied_function(terms, m1_dir, accession):
    """The one function term that the protein's UniProt domains name, if exactly one.
    A kinase domain alone implies 'Kinase activity'; WRN has helicase AND exonuclease
    domains, so nothing is implied there (it stays unlabeled rather than guessed).
    When no domain names a term, the EC number can (EC_IMPLIES): CYP3A4 -> 'Cytochrome P450 activity',
    carbonic anhydrase -> 'Enzyme activity'. Only terms present in the term table are used."""
    try:
        with open(os.path.join(m1_dir, "uniprot_%s.json" % accession), encoding="utf-8") as handle:
            entry = json.load(handle)
    except (OSError, ValueError):
        return None, ""
    features = entry.get("features", [])
    domains = [clean(f.get("description")) for f in features if f.get("type") == "Domain"]
    found = collections.OrderedDict()
    for term in terms:
        if term.axis != "function" or term.applies_to == "cell-based":
            continue
        if term.name in ("Ligand binding affinity", "Covalent adduct formation"):
            continue
        for domain in domains:
            if term.regex.search(norm(domain)):
                found.setdefault(term.name, domain)
    if len(found) == 1:
        name = list(found)[0]
        return name, found[name]
    if not found:                     # no domain names a term: fall back to the EC number
        ecs = ec_numbers(entry)
        name = ec_implied(ecs, {t.name for t in terms if t.axis == "function" and t.applies_to != "cell-based"})
        if name:
            return name, "EC " + ", ".join(ecs)
    return None, ""

# ------------------------------------------------------ language model

class Embedder:
    """Computes sentence vectors with a cache saved on disk."""

    def __init__(self, cache_dir):
        import numpy as np
        from sentence_transformers import SentenceTransformer
        self.np = np
        self.model = SentenceTransformer(MODEL)
        self.keys_path = os.path.join(cache_dir, "embed_keys.json")
        self.vecs_path = os.path.join(cache_dir, "embed_vecs.npy")
        self.index, self.vecs = {}, None
        self.new = 0
        if os.path.exists(self.keys_path) and os.path.exists(self.vecs_path):
            with open(self.keys_path, encoding="utf-8") as handle:
                keys = json.load(handle)
            self.vecs = np.load(self.vecs_path)
            if len(keys) == len(self.vecs):
                self.index = {k: i for i, k in enumerate(keys)}
            else:
                self.vecs = None

    def encode_raw(self, texts):
        return self.model.encode(texts, normalize_embeddings=True, batch_size=64,
                                 show_progress_bar=False)

    def get(self, texts):
        np = self.np
        keys = [hashlib.sha1((MODEL + "|" + t).encode("utf-8")).hexdigest() for t in texts]
        missing = []
        for key, text in zip(keys, texts):
            if key not in self.index and key not in dict(missing):
                missing.append((key, text))
        if missing:
            new_vecs = self.encode_raw([t for _k, t in missing])
            base = 0 if self.vecs is None else len(self.vecs)
            self.vecs = new_vecs if self.vecs is None else np.vstack([self.vecs, new_vecs])
            for i, (key, _t) in enumerate(missing):
                self.index[key] = base + i
            self.new += len(missing)
        return np.array([self.vecs[self.index[k]] for k in keys])

    def save(self):
        if self.vecs is None:
            return
        keys = [None] * len(self.index)
        for key, i in self.index.items():
            keys[i] = key
        with open(self.keys_path, "w", encoding="utf-8") as handle:
            json.dump(keys, handle)
        self.np.save(self.vecs_path, self.vecs)


class SemanticAxis:
    def __init__(self, embedder, terms, axis):
        np = embedder.np
        self.np = np
        self.terms = [t for t in terms if t.axis == axis]
        phrases, starts = [], []
        for term in self.terms:
            starts.append(len(phrases))
            phrases.extend(term.sem_phrases)
        self.starts = np.array(starts)
        self.vecs = embedder.get(phrases) if phrases else None

    def best(self, chunk_vecs, cell_based):
        """Return (best term, score, second term, second score).
        Terms that only apply to cell assays are ruled out for cell-free assays."""
        np = self.np
        if self.vecs is None or not len(self.terms):
            return None, 0.0, None, 0.0
        sims = (chunk_vecs @ self.vecs.T).max(axis=0)
        per_term = np.maximum.reduceat(sims, self.starts).copy()
        if not cell_based:
            for i, term in enumerate(self.terms):
                if term.applies_to == "cell-based":
                    per_term[i] = -1.0
        order = np.argsort(per_term)[::-1]
        first = order[0]
        second = order[1] if len(order) > 1 else None
        return (self.terms[first].name, float(per_term[first]),
                self.terms[second].name if second is not None else None,
                float(per_term[second]) if second is not None else 0.0)


# ------------------------------------------------------------ main work


def label_records(records, terms, embedder, patterns, length, implied=(None, "")):
    sem = {}
    if embedder is not None:
        for axis in ("function", "technology"):
            sem[axis] = SemanticAxis(embedder, terms, axis)
    review = []
    for n, rec in enumerate(records, 1):
        text = rec["text"]
        tnorm = norm(text)
        rec["context"] = context_of(text, rec.get("bao", ""), rec.get("cell_type", ""))
        rule_hits = {axis: rule_match(terms, axis, tnorm) for axis in ("function", "technology")}
        if rec["context"] == "cell-free" and any(
                t.applies_to == "cell-based" and t.name in rule_hits[t.axis] for t in terms):
            rec["context"] = "cell-based"
        if COUNTER.search(text):
            rec["role"] = "counterscreen"
        elif CONFIRM.search(text):
            rec["role"] = "confirmation"
        else:
            rec["role"] = "unspecified"
        rec["effect"] = ("inhibition" if INHIB.search(text) else
                         "activation" if ACTIV.search(text) else "unspecified")
        rec["mutant"] = "yes" if (rec.get("mutation") or MUTANT.search(text)) else ""
        rec["in_v1_scope"] = "yes" if rec["context"] == "cell-free" else "no"
        construct = find_construct(text, patterns, length)
        rec["construct_start"], rec["construct_end"] = construct["start"], construct["end"]
        rec["construct_note"], rec["construct_evidence"] = construct["note"], construct["evidence"]

        chunk_vecs = None
        for axis in ("function", "technology"):
            hits = rule_hits[axis]
            labels = list(hits.keys())
            method = "rule" if labels else ""
            guess = (None, 0.0, None, 0.0)
            if not labels and axis in sem:
                if chunk_vecs is None:
                    chunk_vecs = embedder.get(chunk(text))
                guess = sem[axis].best(chunk_vecs, rec["context"] == "cell-based")
                if guess[0] and guess[1] >= SEM_MIN and guess[1] - guess[3] >= SEM_MARGIN:
                    labels, method = [guess[0]], "semantic %.2f" % guess[1]
            evidence = "; ".join("%s <- '%s'" % kv for kv in hits.items())
            if not labels and axis == "technology" and rec.get("abstract") and rec["context"] == "cell-free":
                # second reading: the paper's abstract, used only when it names exactly one method
                found = rule_match(terms, "technology", norm(rec["abstract"]))
                found = {k: v for k, v in found.items() if k not in ABSTRACT_SKIP
                         and not any(t.name == k and t.applies_to == "cell-based" for t in terms)}
                if len(found) == 1:
                    labels, method = list(found), "abstract"
                    evidence = "; ".join("%s <- '%s' (in paper abstract)" % kv for kv in found.items())
            if (not labels and axis == "function" and implied[0] and rec["context"] == "cell-free"
                    and rec["effect"] != "unspecified" and not BINDING_WORDS.search(text)):
                if implied[1].startswith("EC "):
                    labels, method = [implied[0]], "implied by EC number"
                    evidence = "UniProt %s" % implied[1]
                else:
                    labels, method = [implied[0]], "implied by domain"
                    evidence = "UniProt domain '%s'" % implied[1]
            rec[axis + "_labels"] = " | ".join(labels)
            rec[axis + "_method"] = method
            rec[axis + "_evidence"] = evidence
            if not labels and rec["in_v1_scope"] == "yes":
                review.append({
                    "id": rec["id"], "axis": axis, "text_level": rec["text_level"],
                    "text": clip(text, 260),
                    "best_guess": guess[0] or "", "score": round(guess[1], 2),
                    "second_guess": guess[2] or "", "second_score": round(guess[3], 2),
                    "your_label": "",
                })
        if embedder is not None and n % 200 == 0:
            print("    labeled %d of %d" % (n, len(records)))
    return review


def write_csv(path, rows, columns):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def evenly(rows, n):
    if len(rows) <= n:
        return rows
    step = len(rows) / float(n)
    return [rows[int(i * step)] for i in range(n)]


def summarize(gene, records, dropped, review):
    say("  Merged list: %d assays kept, %d dropped." % (len(records), len(dropped)))
    why = collections.Counter(d["reason"] for d in dropped)
    for reason, count in why.most_common():
        say("    dropped %-34s %d" % (reason, count))
    campaigns = {r["campaign"] for r in records}
    say("  Independent campaigns among kept assays: %d" % len(campaigns))
    say("  Context: " + ", ".join("%s %d" % kv for kv in collections.Counter(
        r["context"] for r in records).most_common()))
    say("  Role guess: " + ", ".join("%s %d" % kv for kv in collections.Counter(
        r["role"] for r in records).most_common()))
    say("  Mutant-form assays: %d" % sum(1 for r in records if r["mutant"]))

    scope = [r for r in records if r["in_v1_scope"] == "yes"]
    say("  In version 1 scope (cell-free): %d assays" % len(scope))
    for axis in ("function", "technology"):
        labeled = sum(1 for r in scope if r[axis + "_labels"])
        how = collections.Counter(
            ("rule" if r[axis + "_method"] == "rule" else
             "meaning" if r[axis + "_method"].startswith("semantic") else r[axis + "_method"])
            for r in scope if r[axis + "_labels"])
        say("  %s labeled: %d of %d (%s); NOT STATED in the text: %d (%d%%)" % (
            axis, labeled, len(scope), ", ".join("%s %d" % kv for kv in how.most_common()),
            len(scope) - labeled, 100 * (len(scope) - labeled) // max(1, len(scope))))

    def count_campaigns(selector):
        found = collections.defaultdict(set)
        for r in scope:
            for label in selector(r):
                found[label].add(r["campaign"])
        return sorted(found.items(), key=lambda kv: -len(kv[1]))

    def per_axis(axis):
        found, by_meaning = collections.defaultdict(set), collections.defaultdict(set)
        for r in scope:
            for label in [x for x in r[axis + "_labels"].split(" | ") if x]:
                found[label].add(r["campaign"])
                if r[axis + "_method"] != "rule":
                    by_meaning[label].add(r["campaign"])
        return sorted(found.items(), key=lambda kv: -len(kv[1])), by_meaning

    for axis, title in (("function", "FUNCTION"), ("technology", "TECHNOLOGY")):
        say("  Campaigns per %s (cell-free); in brackets, how many were inferred (language model, abstract or domain):" % title)
        rows, by_meaning = per_axis(axis)
        for label, camps in rows:
            say("    %3d  (%3d)  %s" % (len(camps), len(by_meaning.get(label, ())), label))

    def pairs(r):
        funcs = [x for x in r["function_labels"].split(" | ") if x] or ["(function unlabeled)"]
        techs = [x for x in r["technology_labels"].split(" | ") if x] or ["(technology unlabeled)"]
        return ["%s + %s" % (f, t) for f in funcs for t in techs]

    say("  Preview ranking: approaches (function + technology) by independent campaigns:")
    for label, camps in count_campaigns(pairs)[:15]:
        say("    %3d  %s" % (len(camps), label))

    ranged = [r for r in records if r["construct_start"] != ""]
    say("  Construct ranges found in text: %d of %d assays" % (len(ranged), len(records)))
    say("  Construct notes: " + ", ".join("%s %d" % kv for kv in collections.Counter(
        r["construct_note"] for r in records).most_common()))
    say("  Left for your review (unlabeled, cell-free): %d rows" % len(review))


def run_gene(gene, args, terms, embedder):
    say("")
    say("=== %s ===" % gene)
    accession = GENES.get(gene.upper())
    if accession is None:
        say("  Unknown gene %s. Known: %s" % (gene, ", ".join(GENES)))
        return
    abstracts = load_abstracts(args.m1)
    records, dropped = build_records(gene, args.m1, abstracts)
    if records is None:
        say("  No M1 files found for %s in %s. Run M1 first." % (gene, args.m1))
        return
    length = protein_length(args.m1, accession)
    say("  Protein length from UniProt: %s" % (length or "unknown"))
    implied = implied_function(terms, args.m1, accession)
    say("  Paper abstracts available: %d papers (%d assays). %s" % (
        len(abstracts), sum(1 for r in records if r.get("abstract")),
        "" if abstracts else "Run m2b_fetch_abstracts.py to add them."))
    say("  Function implied by the UniProt domains or EC number: %s" % (
        "%s (from %s)" % (implied[0], implied[1] if implied[1].startswith("EC ") else "domain '%s'" % implied[1])
        if implied[0] else "none (several or no matching domains)"))
    review = label_records(records, terms, embedder, construct_patterns(gene), length, implied)

    out = args.out
    write_csv(os.path.join(out, "%s_labeled.csv" % gene), records,
              ["id", "source", "campaign", "text_level", "in_v1_scope", "context", "role", "effect",
               "mutant",
               "function_labels", "function_method", "technology_labels", "technology_method",
               "construct_start", "construct_end", "construct_note", "chembl_type", "bao",
               "pubchem_aid", "title", "function_evidence", "technology_evidence",
               "construct_evidence"])
    write_csv(os.path.join(out, "%s_dropped.csv" % gene), dropped,
              ["id", "source", "reason", "title"])
    if len(review) > REVIEW_CAP:
        review = evenly(review, REVIEW_CAP)
    write_csv(os.path.join(out, "%s_review.csv" % gene), review,
              ["id", "axis", "text_level", "text", "best_guess", "score", "second_guess", "second_score",
               "your_label"])
    # spot-check sample: a few rows of every way a label can be made (rule, meaning, abstract, domain)
    spot = []
    for axis in ("function", "technology"):
        for kind in ("rule", "semantic", "abstract", "implied by domain"):
            rows = [r for r in records if r["in_v1_scope"] == "yes" and r[axis + "_labels"]
                    and r[axis + "_method"].startswith(kind)]
            for r in evenly(rows, 8):
                spot.append({"id": r["id"], "axis": axis, "label": r[axis + "_labels"],
                             "made_by": r[axis + "_method"], "text_level": r["text_level"],
                             "evidence": r[axis + "_evidence"], "assay_text": clip(r["text"], 260),
                             "correct_yes_no": ""})
    write_csv(os.path.join(out, "%s_spotcheck.csv" % gene), spot,
              ["id", "axis", "label", "made_by", "text_level", "evidence", "assay_text", "correct_yes_no"])
    checks = evenly([r for r in records if r["construct_start"] != ""], 30)
    write_csv(os.path.join(out, "%s_constructs_check.csv" % gene),
              [{"id": r["id"], "start": r["construct_start"], "end": r["construct_end"],
                "note": r["construct_note"], "text_around_match": r["construct_evidence"],
                "correct_yes_no": ""} for r in checks],
              ["id", "start", "end", "note", "text_around_match", "correct_yes_no"])
    aids = sorted({r["pubchem_aid"] for r in records if r["source"] == "ChEMBL" and r["pubchem_aid"]})
    missing = sorted({r["pubchem_aid"] for r in records if r["source"] == "ChEMBL"
                      and r["pubchem_aid"].isdigit() and r["text_level"] == "ChEMBL one-line description only"})
    with open(os.path.join(out, "%s_original_pubchem_aids.txt" % gene), "w") as handle:
        handle.write("\n".join(missing) + ("\n" if missing else ""))
    summarize(gene, records, dropped, review)
    say("  How much text each assay has: " + ", ".join(
        "%s %d" % kv for kv in collections.Counter(r["text_level"] for r in records).most_common()))
    odd = [a for a in aids if not a.isdigit()]
    if missing:
        say("  Original PubChem deposits still to download (run m1b_fetch_original_aids.py): %s"
            % ", ".join(missing))
    if odd:
        say("  Skipped (not a plain PubChem number): %s" % ", ".join(odd))


def main():
    global SEM_MIN, SEM_MARGIN
    parser = argparse.ArgumentParser()
    parser.add_argument("genes", nargs="*", help="gene symbols, default WRN and EGFR")
    parser.add_argument("--rules-only", action="store_true", help="skip the language model")
    parser.add_argument("--sem-min", type=float, default=None,
                        help="minimum language-model score (default %.2f)" % SEM_MIN)
    parser.add_argument("--sem-margin", type=float, default=None,
                        help="minimum lead over the runner-up (default %.2f)" % SEM_MARGIN)
    parser.add_argument("--reset-terms", action="store_true",
                        help="replace term_table.csv with the newest starter table (old one is kept)")
    parser.add_argument("--m1", default=os.path.join(HERE, "m1_data"))
    parser.add_argument("--out", default=os.path.join(HERE, "m2_data"))
    args = parser.parse_args()
    if args.sem_min is not None:
        SEM_MIN = args.sem_min
    if args.sem_margin is not None:
        SEM_MARGIN = args.sem_margin
    genes = [g.upper() for g in args.genes] or list(GENES)
    os.makedirs(args.out, exist_ok=True)

    say("M2 labeling for the Assay Strategy Planner")
    say("Run at: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    if args.reset_terms and os.path.exists(TERM_TABLE):
        backup = TERM_TABLE.replace(".csv", "_old_%s.csv" % datetime.now().strftime("%Y%m%d_%H%M%S"))
        os.replace(TERM_TABLE, backup)
        say("Old term table kept as: %s" % backup)
    say("Language-model thresholds: score >= %.2f, lead over runner-up >= %.2f" % (SEM_MIN, SEM_MARGIN))
    terms = load_terms()
    say("Term table: %d terms (%d function, %d technology)." % (
        len(terms), sum(t.axis == "function" for t in terms),
        sum(t.axis == "technology" for t in terms)))

    embedder = None
    if not args.rules_only:
        try:
            started = time.time()
            embedder = Embedder(args.out)
            say("Language model ready in %.1f s." % (time.time() - started))
        except Exception as error:
            say("Language model not available (%s: %s). Continuing with synonyms only."
                % (type(error).__name__, error))

    for gene in genes:
        run_gene(gene, args, terms, embedder)

    if embedder is not None:
        embedder.save()
        say("")
        say("Language model cache: %d vectors saved (%d new this run)."
            % (len(embedder.index), embedder.new))

    say("")
    say("DONE. Results are in: %s" % args.out)
    path = os.path.join(HERE, "m2_report.txt")
    try:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(LOG) + "\n")
        say("Report saved as: %s" % path)
    except OSError as error:
        say("Could not save the report (%s). Please copy the text above." % error)


if __name__ == "__main__":
    main()
