#!/usr/bin/env python
"""Is the protein the subject of the paper's assay work, or just the context it was done in?

The problem this solves, observed in the first survey run: MTAP ranked fourth out of 3,441 proteins,
but none of its top papers contains an MTAP assay. They are PRMT5 papers, and MTAP appears only
because MTAP deletion is the patient-stratification marker that makes PRMT5 a target. Mining those
papers for MTAP readouts yields PRMT5 readouts filed under the wrong protein.

The same shape recurs wherever a gene is used as a biomarker rather than assayed:
    "Discovery of selective PRMT5 inhibitors for MTAP-deleted cancers"
    "Synthetic lethality of WRN loss in MSI-high tumours"
    "ATR inhibition in ATM-deficient models"
In each, the gene named second is a condition of the experiment, not its target.

Two decisions are made here:

  mention level   Each occurrence of the protein's name is either a SUBJECT mention or a CONTEXT
                  mention, judged from the words immediately around it. "MTAP-deleted", "loss of
                  MTAP", "MTAP-null", "MTAP status" are context. "MTAP inhibitors", "recombinant
                  MTAP", "the MTAP enzyme" are subject.

  title level     A title is SUBJECT if it has at least one subject mention, CONTEXT if every
                  mention is a context mention, and ABSENT if the name does not appear at all.
                  ABSENT is not a penalty: Europe PMC is queried over title and abstract, so the
                  name may legitimately be in the abstract only.

Two uses:

  survey_papers.py    `subject_fraction` over the returned titles, folded into the supply score, so
                      a protein whose literature is all about somebody else's target ranks down.

  readout_gaps.py     `subject_positions` in place of raw alias positions, so a readout only counts
                      as evidence when it sits near a mention of the protein as an assay subject.

Deliberately conservative. A mention is context only on explicit biomarker wording; anything
ambiguous counts as subject, because the cost of ranking a protein too low is missing real data,
while the cost of ranking it too high is one wasted extraction.

Run this file directly for its self-test, or with --build-index to rebuild the alias index from the
cached UniProt entries (needed once, before the alias guard can judge short abbreviations).
"""

import bisect
import collections
import glob
import json
import os
import re
import sys

# ---------------------------------------------------------------- normalisation

def normalise(text):
    """Same as readout_gaps.normalise and m2_label.normT: lower case, punctuation to single spaces.
    So "MTAP-deleted" becomes "mtap deleted" and one pattern covers hyphen, space and slash."""
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", (text or "").lower()))


def alias_pattern(aliases):
    """Word-bounded alternation over normalised aliases. Boundaries are not optional: without them
    a three-letter gene symbol matches inside ordinary words."""
    alts = sorted({normalise(a) for a in aliases if normalise(a)}, key=len, reverse=True)
    if not alts:
        return None
    return re.compile(r"\b(?:%s)\b" % "|".join(re.escape(a) for a in alts))


# ---------------------------------------------------------------- alias breadth

# Words that name a protein class rather than a protein. An alias built only from these, with no
# number to pin it down, matches every sibling in the family and is useless as a target alias.
FAMILY_WORDS = {
    "receptor", "receptors", "protein", "proteins", "kinase", "kinases", "phosphatase",
    "phosphatases", "protease", "proteases", "peptidase", "transporter", "transporters",
    "channel", "channels", "enzyme", "enzymes", "transferase", "reductase", "oxidase",
    "dehydrogenase", "hydrolase", "esterase", "synthase", "synthetase", "ligase", "isomerase",
    "factor", "factors", "subunit", "subunits", "chain", "domain", "component", "member",
    "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "theta", "like", "type", "family",
    "putative", "probable", "isoform", "isoforms", "variant", "precursor", "mitochondrial",
    "cytosolic", "cytoplasmic", "nuclear", "membrane", "soluble", "secreted", "human",
    "atp", "gtp", "adp", "camp", "cgmp", "nad", "nadp", "purinergic", "serine", "threonine",
    "tyrosine", "dependent", "associated", "binding", "related", "regulated", "activated",
    "and", "of", "the",
}
MIN_ALIAS_LEN = 3   # "TR" and "FP" name nothing; three letters can be real (P53, BRM, PKR, S2R)


# ---- the alias index: how many proteins answer to each name -----------------------------------
#
# Length is the wrong test for whether a short alias is specific. "P53" (TP53), "BRM" (SMARCA2),
# "PKR" (EIF2AK2) and "S2R" (TMEM97) are three letters and each names exactly one protein, while
# "CPA" (CPA1-CPA6), "MST" (MST1, MST1R, MST4), "TIE" (TIE1, TEK) and "CK1" (seven casein kinases)
# are the same length and name a family. A first attempt using a four-character floor dropped P53,
# BRM and PKR, which is worse than the problem it fixed.
#
# So count instead: build an index over the UniProt entries we already have of how many distinct
# proteins list each name. A name claimed by two or more proteins cannot pick out one of them.

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ALIAS_INDEX_PATH = os.path.join(ROOT, "build_cache", "alias_index.json")
_INDEX_CACHE = {}


def entry_aliases(entry):
    """Every name a UniProt entry answers to: gene symbol, gene synonyms, recommended name,
    alternative names. Same extraction readout_gaps.target_aliases uses."""
    out = []
    desc = entry.get("proteinDescription") or {}
    rec = ((desc.get("recommendedName") or {}).get("fullName") or {}).get("value")
    if rec:
        out.append(rec)
    for alt in desc.get("alternativeNames") or []:
        v = (alt.get("fullName") or {}).get("value")
        if v:
            out.append(v)
    for g in entry.get("genes") or []:
        v = ((g.get("geneName") or {}).get("value"))
        if v:
            out.append(v)
        for syn in g.get("synonyms") or []:
            if syn.get("value"):
                out.append(syn["value"])
    return out


def build_alias_index(uniprot_dir, out_path=ALIAS_INDEX_PATH, log=None):
    """{normalised name: number of distinct proteins that answer to it}, over every cached UniProt
    entry. One pass of a few hundred megabytes, then cached as JSON, so this is paid once."""
    counts = collections.Counter()
    symbols = set()
    n = 0
    for path in glob.glob(os.path.join(uniprot_dir, "uniprot_*.json")):
        try:
            with open(path, encoding="utf-8") as h:
                entry = json.load(h)
        except (OSError, ValueError):
            continue
        n += 1
        for name in {normalise(a) for a in entry_aliases(entry)}:
            if name:
                counts[name] += 1
        for g in entry.get("genes") or []:
            sym = normalise((g.get("geneName") or {}).get("value") or "")
            if sym:
                symbols.add(sym)
        if log and n % 1000 == 0:
            log("  indexed %d entries, %d distinct names" % (n, len(counts)))
    # Only ambiguous names are worth storing; a name used once is the common case and the default.
    index = {k: v for k, v in counts.items() if v >= 2}
    if out_path:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as h:
            json.dump({"n_entries": n, "n_names": len(counts), "ambiguous": index,
                       "symbols": sorted(symbols)}, h)
    return index, n


def load_alias_index(path=ALIAS_INDEX_PATH):
    """(ambiguous-name counts, sorted gene symbols), or empty if the index has not been built.
    Absent, alias_is_specific falls back to its wording rules alone, which still catch the
    multi-word family phrases."""
    if path in _INDEX_CACHE:
        return _INDEX_CACHE[path]
    try:
        with open(path, encoding="utf-8") as h:
            d = json.load(h) or {}
        bundle = (d.get("ambiguous") or {}, d.get("symbols") or [])
    except (OSError, ValueError):
        bundle = ({}, [])
    _INDEX_CACHE[path] = bundle
    return bundle


def family_stem_count(norm, symbols):
    """How many gene symbols this name is a strict prefix of. "CPA" is a prefix of CPA1 through
    CPA6, so it is a family stem; "BRM" is a prefix of BRMS1 alone, so it is not. Two or more is
    the threshold, which is what keeps BRM, P53, PKR and S2R."""
    if not symbols or not norm:
        return 0
    i = bisect.bisect_left(symbols, norm)
    n = 0
    while i < len(symbols) and symbols[i].startswith(norm):
        if symbols[i] != norm:
            n += 1
        i += 1
    return n


def alias_is_specific(alias, index=None):
    """Would this SECONDARY alias pick out this protein, or its whole family? The gene symbol itself
    is never put through this test: WRN, ATM, PLG and MET are three letters and all legitimate.

        "TAAR1"                                 yes
        "Trace amine-associated receptor 1"     yes, the number pins it down
        "P53", "BRM", "PKR", "S2R"              yes, one protein each
        "Purinergic receptor"                   NO, matches P2X4 and P2Y2 papers
        "ATP receptor"                          NO
        "CK1"                                   NO, the casein kinase 1 family has seven members
        "CPA", "MST", "TIE"                     NO, family stems
        "TR-beta"                               NO, two letters and a Greek letter
    """
    norm = normalise(alias).strip()
    if len(norm.replace(" ", "")) < MIN_ALIAS_LEN:
        return False
    counts, symbols = load_alias_index() if index is None else index
    if counts.get(norm, 0) >= 2:          # two proteins answer to this name verbatim
        return False
    if family_stem_count(norm, symbols) >= 2:   # it is the stem of a numbered family
        return False
    # A word of one or two letters distinguishes nothing on its own, so it counts as a family word.
    words = [w for w in norm.split() if w]
    if all(w in FAMILY_WORDS or len(w) <= 2 for w in words):
        return False
    if re.search(r"\d", norm):
        return True
    # Without the index, a short all-letters abbreviation is a coin toss, so keep it: a missed
    # sibling paper costs one wasted extraction, a dropped real alias costs the whole protein.
    return True


def specific_aliases(aliases, gene=None, index=None):
    """Keep the aliases that identify the protein; drop the family names.

    `gene` is the gene symbol and is always kept, whatever its length. If it is not given, the
    first alias is assumed to be it. If everything else is dropped the symbol alone remains, so a
    protein never ends up with no name at all."""
    aliases = [a for a in aliases if a]
    if gene is None:
        gene = aliases[0] if aliases else None
    if index is None:
        index = load_alias_index()
    kept = [a for a in aliases if a == gene or alias_is_specific(a, index)]
    if not kept and gene:
        kept = [gene]
    return kept


# ---------------------------------------------------------------- context wording

# Words that, following the protein name, mean the protein is a property of the experimental system
# rather than the thing being measured. Hyphens are already spaces after normalise().
# "dependent" and "dependency" are NOT here: "WRN dependency" means WRN is the target of interest,
# which is the opposite of a biomarker, and excluding them cost real WRN papers in testing.
# "expressing" is NOT here either: "cells expressing TAAR1" is a method describing the assay.
CONTEXT_AFTER = re.compile(
    r"^ ?(?:gene )?(?:"
    r"deleted|deletions?|dels?|delta|deficient|deficiency|null|loss|lost|absent|silenced|silencing|"
    r"hypermethylated|methylated|methylation|knockouts?|ko|kd|negative|positive|low|high|status|"
    r"wild ?type|wt|co ?deleted|co ?deletions?|mutants?|mutated|mutations?|altered|alterations?|"
    r"amplified|amplifications?|proficient|intact|haploinsufficien\w*|nonsense|truncated"
    r")\b")

# The hard genotype words above, used as a document-level signal.
GENOTYPE_HARD = re.compile(
    r"\b(?:wt|wild ?type|dels?|deleted|deletions?|null|knockouts?|ko|deficient|isogenic|"
    r"proficient)\b")

# "<name> cells" / "<name> line". Ambiguous on its own: "HEK293-TAAR1 cells" is the assay system and
# the protein IS the subject, while "HCT116 MTAP cells" is a genotype label on somebody else's model.
# Resolved at document level by `isogenic_pair` below, not here.
CONTEXT_CELLS = re.compile(r"^ ?(?:cells?|cell lines?|lines?|clones?|models?)\b")

# Words that, preceding the protein name as "<word> of <name>", mean the same thing.
CONTEXT_BEFORE = re.compile(
    r"\b(?:"
    r"loss|losses|deletion|deletions|absence|lack|silencing|status|deficiency|methylation|"
    r"depletion|knockdown|knockout|inactivation|downregulation|upregulation|mutations?|"
    r"amplification|expression levels?|levels?"
    r") of $")

# Disease and cell nouns directly after the name. Weaker: "TAAR1 cells" could be an assay system.
# Applied only when another gene in the same title already owns the pharmacological noun, which is
# what makes "PRMT5 inhibitors in MTAP cancers" decidable.
CONTEXT_AFTER_SOFT = re.compile(
    r"^ ?(?:"
    r"cancers?|tumou?rs?|carcinomas?|adenocarcinomas?|melanomas?|gliomas?|glioblastomas?|"
    r"leukemias?|leukaemias?|lymphomas?|sarcomas?|mesotheliomas?|malignanc\w+|neoplas\w+|"
    r"patients?|cohorts?|tissues?|biopsies|xenografts?|models?|backgrounds?|contexts?|settings?"
    r")\b")

# A pharmacological noun bound to a gene-shaped token: "PRMT5 inhibitors", "inhibitors of PRMT5".
PHARM_NOUN = (r"inhibitors?|inhibition|antagonis\w+|agonis\w+|activators?|degraders?|modulators?|"
              r"protacs?|binders?|blockers?|potentiators?|substrates?")
GENE_TOKEN = r"[A-Z][A-Z0-9]{1,}(?:[A-Z0-9-]*[A-Z0-9])?"
OWNS_AFTER = re.compile(r"\b(%s)\b(?:\s*\([^)]{1,20}\))?[\s,/-]{1,3}(?:%s)\b" % (GENE_TOKEN, PHARM_NOUN), re.I)
OWNS_BEFORE = re.compile(r"\b(?:%s)\s+(?:of|against|targeting)\s+(?:the\s+)?(%s)\b" % (PHARM_NOUN, GENE_TOKEN))

# Upper-case tokens that look like gene symbols but are not. Without this stoplist "ATP-competitive
# inhibitors" and "DMSO" make every title look like it is about another protein.
NOT_A_GENE = {
    "DNA", "RNA", "MRNA", "ATP", "ADP", "GTP", "GDP", "AMP", "CAMP", "NAD", "NADH", "NADP", "SAM",
    "FDA", "HTS", "SAR", "NMR", "SPR", "BLI", "ITC", "FRET", "HTRF", "TR", "ELISA", "PCR", "QPCR",
    "LC", "MS", "HPLC", "UPLC", "DMSO", "PDB", "EC", "IC", "KD", "KI", "PK", "PD", "ADME",
    "IN", "VITRO", "VIVO", "SILICO", "AND", "THE", "FOR", "WITH", "VIA", "NEW", "NOVEL", "ONE",
    "II", "III", "IV", "VI", "VII", "ATPASE", "GTPASE", "UV", "IR", "CD", "SPA", "FP", "AI", "ML",
    "SARS", "COV", "HIV", "HCV", "HBV", "TB", "COVID",
}


def _gene_tokens(raw_text, pattern):
    for m in pattern.finditer(raw_text or ""):
        tok = m.group(1)
        if tok.upper() in NOT_A_GENE or not tok.isupper():
            continue
        yield tok


def other_gene_owns_pharmacology(raw_title, aliases):
    """True if a gene-shaped token that is NOT this protein is the one carrying the drug noun.
    Runs on the raw title, not the normalised one, because the capital letters are the signal."""
    mine = {normalise(a) for a in aliases}
    for pattern in (OWNS_AFTER, OWNS_BEFORE):
        for tok in _gene_tokens(raw_title, pattern):
            if normalise(tok) not in mine:
                return True
    return False


# ---------------------------------------------------------------- mention and title classification

def isogenic_pair(norm_text, pat, window=60):
    """Does this document use the protein's name as a cell-line genotype label?

    The giveaway is the isogenic pair: a paper that writes both "HCT116 MTAP WT" and "HCT116 MTAP
    del" is using MTAP status to label two models, not assaying MTAP. Found this in the MTA-
    cooperative PRMT5 papers, where "SDMA ELISA" sits next to "HCT116 MTAP cells" and credited
    ELISA to MTAP. Checked once per document, because it is a property of the document."""
    if pat is None:
        return False
    for m in pat.finditer(norm_text):
        if GENOTYPE_HARD.match(norm_text[m.end():m.end() + window].lstrip()):
            return True
    return False


def mention_is_context(norm_text, start, end, allow_soft=False, genotype_doc=False):
    """Judge one occurrence of the protein's name from the words touching it."""
    after = norm_text[end:end + 34]
    if CONTEXT_AFTER.match(after):
        return True
    if allow_soft and CONTEXT_AFTER_SOFT.match(after):
        return True
    if genotype_doc and CONTEXT_CELLS.match(after):
        return True
    return bool(CONTEXT_BEFORE.search(norm_text[max(0, start - 40):start]))


def subject_positions(norm_text, pat, allow_soft=False, genotype_doc=None):
    """Start offsets of the mentions that read as assay subject. Drop-in replacement for the raw
    alias positions used by a proximity test.

    `genotype_doc` defaults to detecting the isogenic pair in this text. Pass False to disable."""
    if pat is None:
        return []
    if genotype_doc is None:
        genotype_doc = isogenic_pair(norm_text, pat)
    return [m.start() for m in pat.finditer(norm_text)
            if not mention_is_context(norm_text, m.start(), m.end(), allow_soft, genotype_doc)]


def classify_title(title, aliases, pat=None, gene=None):
    """"subject", "context" or "absent" for one title. Family-name aliases are dropped first, so a
    paper about a sibling protein does not read as a paper about this one."""
    aliases = specific_aliases(list(aliases), gene)
    pat = pat or alias_pattern(aliases)
    if pat is None:
        return "absent"
    norm = normalise(title)
    hits = list(pat.finditer(norm))
    if not hits:
        return "absent"
    soft = other_gene_owns_pharmacology(title, aliases)
    for m in hits:
        if not mention_is_context(norm, m.start(), m.end(), allow_soft=soft):
            return "subject"
    return "context"


def subject_fraction(titles, aliases, prior=2.0, gene=None):
    """How much of this protein's literature is actually about this protein.

    Returns (fraction, n_subject, n_context, n_absent). The fraction is shrunk toward 1.0 by
    `prior` pseudo-counts, so one biomarker title out of two does not halve a protein's rank
    while six out of six still does most of the damage:

        0 of 0 observed   -> 1.00   nothing learned, no penalty
        6 subject, 0 ctx  -> 1.00
        3 subject, 3 ctx  -> 0.63
        0 subject, 6 ctx  -> 0.25
    """
    n_sub = n_ctx = n_abs = 0
    aliases = specific_aliases(list(aliases), gene)
    pat = alias_pattern(aliases)
    for t in titles or []:
        verdict = classify_title(t.get("title") if isinstance(t, dict) else t, aliases, pat, gene)
        if verdict == "subject":
            n_sub += 1
        elif verdict == "context":
            n_ctx += 1
        else:
            n_abs += 1
    frac = (n_sub + prior) / (n_sub + n_ctx + prior) if (n_sub + n_ctx + prior) else 1.0
    return round(frac, 3), n_sub, n_ctx, n_abs


# ---------------------------------------------------------------- self-test

CASES = [
    # (title, gene, extra aliases, expected)
    ("Discovery of selective PRMT5 inhibitors for the treatment of MTAP-deleted cancers",
     "MTAP", [], "context"),
    ("MTAP deletion confers enhanced dependency on PRMT5 in cancer cells", "MTAP", [], "context"),
    ("Loss of MTAP sensitises tumours to PRMT5 inhibition", "MTAP", [], "context"),
    ("Targeting PRMT5 in MTAP-null glioblastoma", "MTAP", [], "context"),
    ("A continuous spectrophotometric assay for MTAP phosphorylase activity",
     "MTAP", ["methylthioadenosine phosphorylase"], "subject"),
    ("Transition state analogue inhibitors of human MTAP", "MTAP", [], "subject"),
    # WRN: the real papers must survive
    ("Fragment-based discovery of WRN helicase inhibitors", "WRN", ["werner syndrome atp dependent helicase"], "subject"),
    ("An ADP-Glo assay for WRN ATPase activity in high-throughput screening", "WRN", [], "subject"),
    ("Synthetic lethality of WRN loss in microsatellite-unstable cancers", "WRN", [], "context"),
    ("CRISPR screening identifies modulators of WRN dependency", "WRN", [], "subject"),
    # ATM/ATR: both genes real targets, one used as the context
    ("ATR inhibition is selectively lethal in ATM-deficient tumour cells", "ATM", [], "context"),
    ("ATR inhibition is selectively lethal in ATM-deficient tumour cells", "ATR", [], "subject"),
    # ordinary assay titles
    ("Development of an HTRF assay for EGFR kinase inhibitors", "EGFR", [], "subject"),
    ("A GloSensor cAMP assay for TAAR1 agonists", "TAAR1", [], "subject"),
    ("Identification of positive allosteric modulators of GLP-1R", "GLP1R", ["glp 1r"], "subject"),
    ("Caspase-4 activation in the non-canonical inflammasome pathway", "CASP4", ["caspase 4"], "subject"),
    ("Prognostic value of CASP4 expression in colorectal cancer patients", "CASP4", [], "subject"),
    ("EGFR status predicts response to PRMT5 inhibitors", "EGFR", [], "context"),
    ("Structure-based design of BCHE reactivators", "BCHE", [], "subject"),
    ("Pyrazole inhibitors of ATP-competitive kinases", "MTAP", [], "absent"),
    # Over-broad aliases, from the real UniProt synonym lists. Without the alias guard each of
    # these reads as a paper about the target.
    ("Purinergic Receptor P2Y2 Regulates Neurogenesis Following Spinal Cord Injury",
     "P2RX3", ["P2X purinoceptor 3", "ATP receptor", "Purinergic receptor"], "absent"),
    ("New hemiaminal ether dihydropyrimidinones as CK1/CK2 inhibitors targeting alpha-synuclein",
     "CSNK1A1L", ["Casein kinase I isoform alpha-like", "CK1"], "absent"),
    # ...while the specific aliases still work
    ("A radioligand binding assay for the P2X purinoceptor 3 receptor",
     "P2RX3", ["P2X purinoceptor 3", "ATP receptor", "Purinergic receptor"], "subject"),
]

# (alias, is it specific to one protein?) Needs the alias index: run --build-index first.
ALIAS_CASES = [
    # full names, unambiguous
    ("TAAR1", True), ("Trace amine-associated receptor 1", True),
    ("Thioredoxin reductase 2, mitochondrial", True), ("Selenoprotein Z", True),
    ("Cyclin-dependent kinase-like 5", True), ("P2X purinoceptor 3", True),
    # short abbreviations that name exactly one protein, and must survive
    ("P53", True), ("BRM", True), ("PKR", True), ("S2R", True),
    # family phrases, caught by wording
    ("Purinergic receptor", False), ("ATP receptor", False), ("protein kinase", False),
    ("serine/threonine-protein kinase", False), ("Peptidase M", False), ("TR-beta", False),
    # claimed verbatim by more than one protein. "MSAP" is listed under MTAP, and the index found a
    # second protein in our own set claiming it, which is exactly the case length cannot detect.
    ("CK1", False), ("MSAP", False),
    # family stems: prefix of two or more gene symbols
    ("CPA", False), ("MST", False), ("MGA", False),
]


def _self_test():
    bad = 0
    print("alias specificity")
    for alias, want in ALIAS_CASES:
        got = alias_is_specific(alias)
        if got != want:
            bad += 1
        print("  %s %-40s want %-5s got %s"
              % ("ok  " if got == want else "FAIL", alias, want, got))

    print("\ntitle classification")
    for title, gene, extra, want in CASES:
        got = classify_title(title, [gene] + extra, gene=gene)
        flag = "ok  " if got == want else "FAIL"
        if got != want:
            bad += 1
        print("  %s %-9s want %-8s got %-8s %s" % (flag, gene, want, got, title[:66]))
    print("\n%d case(s) wrong out of %d" % (bad, len(CASES) + len(ALIAS_CASES)))

    print("\ncell-line genotype labels in a methods section")
    prmt5 = ("cells hct116 mtap wt and hct116 mtap del were grown in dmem tumors from both sides "
             "were collected for sdma elisa and hct116 mtap cells were implanted bilaterally")
    taar1 = ("hek293 taar1 cells were plated at 20 000 per well and camp accumulation in taar1 "
             "cells was measured with the glosensor assay")
    for name, text, gene in (("PRMT5 paper, MTAP", prmt5, "mtap"), ("TAAR1 paper, TAAR1", taar1, "taar1")):
        pat = alias_pattern([gene])
        iso = isogenic_pair(text, pat)
        n_all = len(pat.findall(text))
        n_sub = len(subject_positions(text, pat))
        verdict = "ok  " if (iso and n_sub == 0) or (not iso and n_sub == n_all) else "FAIL"
        if verdict == "FAIL":
            bad += 1
        print("  %s %-20s isogenic_pair=%-5s mentions %d, subject %d"
              % (verdict, name, iso, n_all, n_sub))

    print("\nsubject_fraction examples")
    for name, titles, gene in [
        ("MTAP, all PRMT5 papers", [c[0] for c in CASES[:4]], "MTAP"),
        ("MTAP, mixed", [c[0] for c in CASES[:6]], "MTAP"),
        ("WRN", [c[0] for c in CASES[6:10]], "WRN"),
    ]:
        print("   %-26s %s" % (name, subject_fraction(titles, [gene], gene=gene)))
    return 1 if bad else 0


def _build():
    d = os.path.join(ROOT, "build_cache", "m1_data")
    print("Indexing UniProt names in %s ..." % os.path.relpath(d, ROOT))
    index, n = build_alias_index(d, log=lambda m: print(m, flush=True))
    print("  %d entries, %d ambiguous name(s) written to %s"
          % (n, len(index), os.path.relpath(ALIAS_INDEX_PATH, ROOT)))
    _INDEX_CACHE.clear()
    counts, symbols = load_alias_index()
    print("  %-22s %8s %8s %s" % ("name", "claimed", "stem of", "kept?"))
    for name in ("p53", "brm", "pkr", "s2r", "cpa", "mst", "tie", "ck1", "ta1", "mga",
                 "purinergic receptor", "atp receptor", "peptidase m"):
        print("  %-22s %8d %8d %s"
              % (name, counts.get(name, 1), family_stem_count(name, symbols),
                 "keep" if alias_is_specific(name) else "DROP"))
    return 0


if __name__ == "__main__":
    if "--build-index" in sys.argv:
        raise SystemExit(_build())
    raise SystemExit(_self_test())
