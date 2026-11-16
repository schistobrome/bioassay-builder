#!/usr/bin/env python
"""Find, fetch, score and split open-access papers that describe assays for a protein, then (optionally)
fill schema/paper_extraction.schema.json for each assay block, with a local LLM or with the agent.

Europe PMC is used rather than PubChem or PubMed: it is the only one of the three whose API returns
open-access FULL TEXT, and its section-scoped fields (TITLE_ABS:, METHODS:) are what make the queries
selective. PubChem holds no full text at all.

Four stages, each cached on disk so nothing is downloaded twice:

  search   Europe PMC query per protein. Target names go in TITLE_ABS so the paper is ABOUT the protein;
           assay terms go in METHODS so the match is a real protocol, not a passing mention. Terms come
           from data/paper_queries.json and cover inhibitors AND activators (agonist / PAM / potentiator /
           EC50 / Emax), because 11 of the 73 readable sample papers were activation studies.
  fetch    Full text XML for each hit with a PMCID, into build_cache/papers/fulltext/.
  score    Regex scoring for assay-detail density, because relevance order alone is not enough: for CDK2
           it put a cancer-biology paper above the allosteric-inhibitor assay paper. Signals and weights
           live in data/paper_queries.json. Checked against the 75 sample PDFs: 22 scored 5-6 and all 22
           were small-molecule assay papers, mean 3.6 vs 2.0 for biologics.
  split    Cut the methods into one block per subsection heading (<sec><title>), which is the natural
           assay boundary, and write build_cache/papers/blocks/<pmcid>/<n>_<slug>.txt.

Then --extract fills the JSON for each block, by one of three routes:

  --backend ollama   a local model over /api/generate (default)
  --backend openai   any OpenAI-compatible server: llama.cpp --server, vLLM, LM Studio
  --backend claude   write a work queue for the Claude Code agent in THIS session to fill by hand,
                     then assemble it with --agent-merge. The script cannot call Claude Code, because
                     Claude Code is what runs the script, so this inverts the direction: the script
                     prepares and validates, the agent reads the blocks and writes the JSON. Slower
                     and costlier per paper than a local model, and better, so it suits a handful of
                     high-value papers rather than a sweep over thousands.

Whichever route, the finished record is validated against schema/paper_extraction.schema.json.

Examples
--------
  # what would be searched, no network
  py scripts\\fetch_papers.py --genes WRN,EGFR --dry-run

  # search, fetch, score, split
  py scripts\\fetch_papers.py --genes WRN,EGFR,GLP1R --keep 10

  # every protein the app knows, 10 papers each, resumable
  py scripts\\fetch_papers.py --all --keep 10 --max-minutes 120

  # then extract with a local model
  py scripts\\fetch_papers.py --genes WRN --extract --model qwen2.5:7b-instruct

  # or hand the blocks to the Claude Code agent, fill them, merge and validate
  py scripts\\fetch_papers.py --genes WRN --keep 2 --extract --backend claude
  #   ... agent works build_output/papers/agent_queue/TASKS.md ...
  py scripts\\fetch_papers.py --agent-merge
  py scripts\\fetch_papers.py --validate build_output/papers

Standard library only for stages 1-4. --extract needs `jsonschema` to validate, and a local LLM server.
"""

import argparse
import collections
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
SCHEMA = os.path.join(ROOT, "schema", "paper_extraction.schema.json")
CACHE = os.path.join(ROOT, "build_cache", "papers")
SEARCH_DIR = os.path.join(CACHE, "search")
FULL_DIR = os.path.join(CACHE, "fulltext")
BLOCK_DIR = os.path.join(CACHE, "blocks")
OUT_DIR = os.path.join(ROOT, "build_output", "papers")

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
UA = "bioassay-builder literature miner (research use)"
PAUSE = [1.0]

sys.path.insert(0, HERE)
import target_subject as tsub        # noqa: E402

# A title whose only mention of the protein is a biomarker mention, such as "...for MTAP-deleted
# cancers". The paper is not discarded on the title alone, because the methods can still contain a
# real assay for the protein, but it must out-score the papers that are actually about it.
CONTEXT_TITLE_PENALTY = -3

# Methods-like section titles. Everything else in the body is skipped.
METHODS_RE = re.compile(r"method|experimental|materials|protocol", re.I)
# Subsections that are never an assay, so the model is not spent on them.
SKIP_BLOCK_RE = re.compile(
    r"^(statistic|data availability|code availability|reporting summary|ethic|study approval|"
    r"patient|animal|mice|mouse husbandry|xenograft|in vivo|pharmacokinetic|synthesis|general chemistry|"
    r"preparation of|crystal (growth|data)|structure determination|refinement|data collection|"
    r"molecular (docking|dynamics)|ligand docking|sequence alignment|phylogen|rna-?seq|"
    r"next generation sequencing|accession|acknowledg|author contribution|supplementary)", re.I)
ASSAY_HINT_RE = re.compile(
    r"assay|screen|binding|activity|titration|resonance|interferometry|calorimetry|unwinding|"
    r"viability|proliferation|cytotox|imaging|cytometry|blot|elisa|expression and purification|"
    r"protein (expression|production|purification)|reporter|camp|calcium|patch|fret|polarisation|"
    r"polarization|thermal|engagement|degradation|potency|inhibition|potentiat|agonis", re.I)


def log(msg):
    print(msg)
    sys.stdout.flush()


def read_json(path, default=None):
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return default


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=1, ensure_ascii=False)


def http_get(url, timeout=90, tries=3):
    last = None
    for n in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (400, 404):
                raise
            time.sleep(2 ** n)
    raise last


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))[:80]


# ------------------------------------------------------------------ config and target names

def load_config():
    cfg = read_json(os.path.join(DATA, "paper_queries.json"))
    if not cfg:
        sys.exit("Cannot read data/paper_queries.json")
    return cfg


def read_csv_rows(path):
    import csv
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


# Synonyms that are ordinary English words or too generic; they would match thousands of unrelated papers.
STOP_SYNONYMS = {
    "max", "cat", "set", "rest", "impact", "spin", "clock", "hairless", "lamp", "mark", "pigs", "sds",
    "arc", "ash", "bad", "bid", "bod", "cast", "cps", "dad", "egg", "fat", "gap", "hat", "its", "jun",
    "key", "lag", "map", "mat", "mel", "mia", "net", "nut", "pet", "pip", "ran", "sat", "sel", "ski",
    "son", "tag", "tap", "tip", "tra", "van", "was", "win", "ace", "age", "aim", "air", "all", "and",
}
# Clinical aliases that methods sections actually use but UniProt may not list as a gene synonym.
EXTRA_ALIASES = {
    "ERBB2": ["HER2"], "KCNH2": ["hERG"], "GLP1R": ["GLP-1R"], "EGFR": ["HER1"],
    "WRN": ["Werner syndrome helicase"], "PDCD1": ["PD-1"], "CD274": ["PD-L1"],
    "HSP90AA1": ["HSP90"], "BIRC5": ["survivin"], "FOLH1": ["PSMA"],
}


def uniprot_names(accession, gene):
    """Gene symbol plus UniProt names, from the entry the pipeline already cached. No network."""
    names = []
    if gene:
        names.append(gene.upper())
    entry = read_json(os.path.join(ROOT, "build_cache", "m1_data", "uniprot_%s.json" % accession), {}) or {}
    desc = entry.get("proteinDescription") or {}
    rec = ((desc.get("recommendedName") or {}).get("fullName") or {}).get("value")
    if rec:
        names.append(rec)
    for alt in (desc.get("alternativeNames") or [])[:2]:
        v = (alt.get("fullName") or {}).get("value")
        if v:
            names.append(v)
    for g in entry.get("genes") or []:
        for syn in (g.get("synonyms") or [])[:3]:
            v = (syn.get("value") or "").strip()
            if v and len(v) >= 3 and v.lower() not in STOP_SYNONYMS:
                names.append(v)
    names += EXTRA_ALIASES.get((gene or "").upper(), [])
    out, seen = [], set()
    for n in names:
        n = " ".join(str(n).split())
        if n and len(n) >= 3 and n.lower() not in seen and n.lower() not in STOP_SYNONYMS:
            seen.add(n.lower())
            out.append(n)
    return out[:6]


def classes_for(accession, by_acc):
    return by_acc.get(accession, [])


def build_query(cfg, names, class_ids):
    terms = list(cfg["assay_terms_core"])
    for c in class_ids:
        terms += cfg["assay_terms_by_class"].get(c, [])
    seen, uniq = set(), []
    for t in terms:
        if t not in seen:
            seen.add(t)
            uniq.append(t)
    name_field = cfg.get("target_name_field", "TITLE_ABS")
    term_field = cfg.get("assay_term_field", "METHODS")
    name_q = " OR ".join('%s:"%s"' % (name_field, n) for n in names)
    term_q = " OR ".join("%s:%s" % (term_field, t) for t in uniq)
    return cfg["query_template"].format(target_names=name_q, assay_terms=term_q,
                                        exclusions=cfg.get("exclusions", "")).strip()


# ------------------------------------------------------------------ stage 1: search

def search(query, page_size, cfg):
    params = {"query": query, "format": "json", "resultType": "core", "pageSize": min(page_size, 100)}
    if cfg.get("sort"):
        params["sort"] = cfg["sort"]
    data = json.loads(http_get(EPMC + "/search?" + urllib.parse.urlencode(params)))
    return data.get("hitCount", 0), (data.get("resultList") or {}).get("result") or []


# ------------------------------------------------------------------ stage 2: fetch full text

def fetch_fulltext(pmcid):
    path = os.path.join(FULL_DIR, pmcid + ".xml")
    if os.path.exists(path) and os.path.getsize(path) > 500:
        with open(path, "rb") as handle:
            return handle.read()
    xml = http_get("%s/%s/fullTextXML" % (EPMC, pmcid))
    os.makedirs(FULL_DIR, exist_ok=True)
    with open(path, "wb") as handle:
        handle.write(xml)
    time.sleep(PAUSE[0])
    return xml


def el_text(el):
    """itertext(), not .text: on one test paper .text gave 24,754 characters and itertext() gave 78,240,
    because the text after nested italic / sup tags is in their tails."""
    return " ".join("".join(el.itertext()).split())


def body_text(xml):
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return "", None
    body = root.find(".//body")
    return (el_text(body) if body is not None else ""), root


# ------------------------------------------------------------------ stage 3: score

def compile_signals(cfg):
    sig = []
    for name, spec in cfg["ranking"]["score_signals"].items():
        sig.append((name, re.compile(spec["pattern"], re.I), spec.get("weight", 1)))
    pen = cfg["ranking"].get("penalties", {}).get("biologic_only")
    pen_re = re.compile(pen["pattern"], re.I) if pen else None
    pen_w = pen.get("weight", -3) if pen else 0
    small_re = re.compile(r"small[- ]molecule|compound (library|collection)|\bic50\b|\bec50\b|dmso stock|"
                          r"serial dilution|fragment (screen|library)|\bsar\b", re.I)
    return sig, pen_re, pen_w, small_re


def score_text(text, sig, pen_re, pen_w, small_re):
    hits, total = {}, 0
    for name, rx, weight in sig:
        n = len(rx.findall(text))
        hits[name] = n
        if n:
            total += weight
    if pen_re is not None and pen_re.search(text) and not small_re.search(text):
        total += pen_w
        hits["biologic_only_penalty"] = 1
    return total, hits


# ------------------------------------------------------------------ stage 4: split into assay blocks

def assay_blocks(root):
    """One block per methods subsection. Falls back to the whole methods section when it has no subsections."""
    if root is None:
        return []
    body = root.find(".//body")
    if body is None:
        return []
    out = []
    for sec in body.iter("sec"):
        title_el = sec.find("title")
        title = el_text(title_el) if title_el is not None else ""
        if not METHODS_RE.search(title):
            continue
        subs = sec.findall("sec")
        if not subs:
            whole = el_text(sec)
            if len(whole) > 400:
                out.append((title or "Methods", strip_title(whole, title)))
            continue
        for sub in subs:
            st = sub.find("title")
            name = el_text(st) if st is not None else "(untitled)"
            text = strip_title(el_text(sub), name)
            if len(text) < 200:
                continue
            if SKIP_BLOCK_RE.search(name) or not ASSAY_HINT_RE.search(name + " " + text[:400]):
                continue
            out.append((name, text))
    seen, uniq = set(), []
    for name, text in out:
        key = (name.lower(), len(text))
        if key not in seen:
            seen.add(key)
            uniq.append((name, text))
    return uniq


def strip_title(text, title):
    """A subsection body repeats its own title: 'ATPase assayTo measure WRN ATPase activity...'."""
    t = title.strip()
    if t and text.startswith(t):
        return text[len(t):].strip()
    return text


# ------------------------------------------------------------------ stage 5: local LLM extraction

PROMPT_HEADER = """You extract assay methods from biomedical papers into strict JSON.

Return ONE JSON object for the ONE assay described in the METHODS BLOCK below. No prose, no markdown
fence, JSON only.

Hard rules:
1. Never guess. If the block does not state something, use null, and add that field's dotted path to
   "not_reported". A null you cannot justify is worse than a missing field.
2. Concentrations must be the FINAL in-assay value, not the stock. "compounds dissolved in 100% DMSO"
   is a stock, so conditions.dmso_pct stays null; "1% DMSO final" gives 1.
3. Only use a readout, buffer or concentration that belongs to THIS assay. Papers describe several
   assays; do not borrow a neighbouring assay's buffer or a binding panel's readout.
4. A residue range needs an explicit cue ("residues 500-942", "helicase core (500-942)"). A bare
   "(2015-2018)" or "(12-17)" is a citation year or page range: leave construct_start null.
5. If the block says the method was "as described previously", add "method_by_reference" to
   extraction_flags and keep only what THIS block states.
6. Cover activators as carefully as inhibitors. If the assay reads activation, set effect_direction to
   "activation" (or "both"), fill "pharmacology" (agonist / positive_allosteric_modulator / ...), and
   fill "stimulation" with the reference agonist and its level (EC20, EC80, ECmax, unstimulated).
   An EC50 from a potentiator screen is meaningless without the agonist challenge level.
7. Put every term you cannot map to the allowed values in "unmapped_terms" rather than forcing it into
   the nearest one.
8. "evidence" must hold at least one EXACT quote copied from the block, with the field it supports.

Allowed function_labels (choose only from these, or leave the array empty):
__FUNCTION_LABELS__

Allowed technology_labels (choose only from these, or leave the array empty):
__TECHNOLOGY_LABELS__

Allowed assay_category: biochemical_hts, biophysical_binding, cell_functional,
mechanistic_counterscreen, orthogonal_confirmation
Allowed scope: cell-free, cell-based
Allowed role: primary, confirmation, counterscreen, orthogonal, selectivity_panel, unclear
Allowed effect_direction: inhibition, activation, both, binding_only, unclear
Allowed pharmacology: agonist, partial_agonist, inverse_agonist, antagonist,
positive_allosteric_modulator, negative_allosteric_modulator, allosteric_agonist,
orthosteric_inhibitor, allosteric_inhibitor, covalent_inhibitor, degrader, stabiliser,
direct_enzyme_activator, disruptor, unclear
Allowed stimulation.agonist_level: unstimulated/basal, EC20, EC50, EC80, ECmax, saturating, other

Fixed-unit condition keys (give a bare number in that unit, or null):
ph, incubation_min, dmso_pct, mgcl2_mm, nacl_mm, atp_um, bsa_pct, volume_ul, buffer_mm,
detergent_pct, reducing_agent_mm, kcl_mm, mncl2_mm, edta_egta_mm, glycerol_pct, preincubation_min
Text condition keys: temperature, buffer, detergent, reducing_agent, carrier_protein, quench

Shape (omit nothing at the top level; use null or [] where there is nothing):
{"assay_id": "__ASSAY_ID__", "source_section": "__SECTION__", "assay_name_verbatim": "...",
 "assay_category": null, "scope": "cell-free", "role": null,
 "function_labels": [], "technology_labels": [], "effect_direction": null,
 "pharmacology": [], "throughput": null,
 "stimulation": null, "protocol_additions": [],
 "target": null, "substrate": null, "conditions": null, "plate": null, "cells": null,
 "instrumentation": null, "compound_handling": null, "quality": null,
 "screening_campaign": null, "controls": null, "interference_and_counterscreens": null,
 "outputs": [], "analysis": null,
 "evidence": [{"field": "...", "quote": "...", "section": "__SECTION__"}],
 "not_reported": [], "confidence": "high", "extraction_flags": [],
 "unmapped_terms": []}

METHODS BLOCK (paper __PMCID__, section "__SECTION__"):
---
__BLOCK_TEXT__
---
JSON:"""


def term_lists():
    rows = read_csv_rows(os.path.join(HERE, "term_table.csv"))
    fn = [r["term"].strip() for r in rows if (r.get("axis") or "").strip().lower() == "function" and r.get("term")]
    tech = [r["term"].strip() for r in rows if (r.get("axis") or "").strip().lower() == "technology" and r.get("term")]
    return fn, tech


def build_prompt(pmcid, section, text, assay_id, fn, tech, max_chars):
    """Sentinel replacement, not %-formatting or .format(): the prompt contains literal '%' signs
    ("100% DMSO") and literal JSON braces, so both of those would raise."""
    out = PROMPT_HEADER
    for token, value in (("__FUNCTION_LABELS__", "; ".join(fn)),
                         ("__TECHNOLOGY_LABELS__", "; ".join(tech)),
                         ("__ASSAY_ID__", assay_id),
                         ("__PMCID__", pmcid),
                         ("__SECTION__", section),
                         ("__BLOCK_TEXT__", text[:max_chars])):
        out = out.replace(token, value)
    return out


def call_ollama(host, model, prompt, timeout, num_ctx):
    body = json.dumps({
        "model": model, "prompt": prompt, "stream": False, "format": "json",
        "options": {"temperature": 0, "num_ctx": num_ctx},
    }).encode("utf-8")
    req = urllib.request.Request(host.rstrip("/") + "/api/generate", data=body,
                                 headers={"Content-Type": "application/json", "User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read()).get("response", "")


def call_openai_compatible(host, model, prompt, timeout, api_key=None):
    body = json.dumps({
        "model": model, "temperature": 0,
        "response_format": {"type": "json_object"},
        "messages": [{"role": "user", "content": prompt}],
    }).encode("utf-8")
    headers = {"Content-Type": "application/json", "User-Agent": UA}
    if api_key:
        headers["Authorization"] = "Bearer " + api_key
    req = urllib.request.Request(host.rstrip("/") + "/v1/chat/completions", data=body, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        d = json.loads(resp.read())
    return d["choices"][0]["message"]["content"]


def parse_model_json(raw):
    raw = (raw or "").strip()
    if raw.startswith("```"):
        raw = re.sub(r"^```[a-zA-Z]*\s*", "", raw)
        raw = re.sub(r"```\s*$", "", raw).strip()
    try:
        return json.loads(raw), None
    except ValueError as exc:
        m = re.search(r"\{.*\}", raw, re.S)
        if m:
            try:
                return json.loads(m.group(0)), None
            except ValueError as exc2:
                return None, str(exc2)
        return None, str(exc)


# ------------------------------------------------------------------ driver

def protein_rows(args):
    """(accession, gene) for the proteins asked for, and the classes the app assigned them."""
    rows = read_csv_rows(os.path.join(ROOT, "build_cache", "uniprot_human_reviewed.csv"))
    by_gene = {}
    for r in rows:
        g = (r.get("gene") or "").strip().upper()
        if g and g not in by_gene:
            by_gene[g] = r.get("accession")
    by_acc = {}
    summary_dir = os.path.join(ROOT, "build_cache", "summary")
    if os.path.isdir(summary_dir):
        for f in os.listdir(summary_dir):
            if f.endswith(".json"):
                s = read_json(os.path.join(summary_dir, f)) or {}
                if s.get("acc"):
                    by_acc[s["acc"]] = s.get("classes") or []
    wanted = []
    if args.genes:
        for g in [x.strip().upper() for x in args.genes.split(",") if x.strip()]:
            wanted.append((by_gene.get(g, g), g))
    elif args.all:
        seen = set()
        for f in sorted(os.listdir(summary_dir) if os.path.isdir(summary_dir) else []):
            if not f.endswith(".json"):
                continue
            s = read_json(os.path.join(summary_dir, f)) or {}
            g = (s.get("gene") or "").upper()
            if g and g not in seen:
                seen.add(g)
                wanted.append((s.get("acc"), g))
    return wanted, by_acc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", help="comma-separated gene symbols, e.g. WRN,EGFR,GLP1R")
    ap.add_argument("--all", action="store_true", help="every protein with a cached summary")
    ap.add_argument("--pool", type=int, default=None, help="candidates to consider per protein (default from config: 50)")
    ap.add_argument("--keep", type=int, default=None, help="papers to keep per protein (default from config: 10)")
    ap.add_argument("--min-score", type=int, default=None, help="minimum assay-detail score to split and extract")
    ap.add_argument("--dry-run", action="store_true", help="print the queries and stop; no network")
    ap.add_argument("--refresh", action="store_true", help="re-run searches even when a cached result exists")
    ap.add_argument("--max-minutes", type=float, default=0, help="stop after this long (resumable)")
    ap.add_argument("--pause", type=float, default=1.0, help="seconds between API calls (be polite)")
    ap.add_argument("--extract", action="store_true", help="fill the JSON for each assay block")
    ap.add_argument("--backend", choices=["ollama", "openai", "claude"], default="ollama",
                    help="ollama (default, /api/generate); openai (any OpenAI-compatible /v1/chat/completions); "
                         "claude = write a work queue for the Claude Code agent in this session to fill by hand, "
                         "then merge it with --agent-merge")
    ap.add_argument("--agent-merge", action="store_true",
                    help="collect the per-assay JSON files written by the agent into one record per paper, "
                         "validate, and write build_output/papers/<pmcid>.json")
    ap.add_argument("--validate", metavar="PATH", default=None,
                    help="validate one record, or every *.json in a directory, against the schema and stop")
    ap.add_argument("--host", default=None, help="default http://localhost:11434 for ollama, http://localhost:8000 for openai")
    ap.add_argument("--model", default="qwen2.5:14b-instruct", help="local model name")
    ap.add_argument("--api-key", default=os.environ.get("LOCAL_LLM_API_KEY"), help="only if your local server wants one")
    ap.add_argument("--num-ctx", type=int, default=16384, help="context window for ollama")
    ap.add_argument("--block-chars", type=int, default=12000, help="characters of a block sent to the model")
    ap.add_argument("--llm-timeout", type=float, default=600)
    ap.add_argument("--max-blocks", type=int, default=0, help="cap blocks per paper (0 = all)")
    args = ap.parse_args(argv)

    if args.validate:
        return validate_path(args.validate)
    if args.agent_merge:
        return agent_merge()
    if not args.genes and not args.all:
        ap.error("give --genes WRN,EGFR or --all, or --agent-merge, or --validate PATH")
    PAUSE[0] = args.pause
    cfg = load_config()
    pool = args.pool or cfg["ranking"]["candidate_pool"]
    keep = args.keep or cfg["ranking"]["keep"]
    min_score = args.min_score if args.min_score is not None else cfg["ranking"]["minimum_score_to_extract"]
    deadline = time.time() + args.max_minutes * 60 if args.max_minutes else None

    wanted, by_acc = protein_rows(args)
    if not wanted:
        sys.exit("No proteins resolved. Use --genes, or run build_database.py first so summaries exist.")
    log("fetch_papers.py: %d protein(s); pool %d, keep %d, min score %d" % (len(wanted), pool, keep, min_score))

    sig, pen_re, pen_w, small_re = compile_signals(cfg)
    fn_terms, tech_terms = term_lists()
    if args.extract and not fn_terms:
        log("  WARNING: scripts/term_table.csv not readable, so the prompt has no vocabulary list.")
    host = args.host or ("http://localhost:11434" if args.backend == "ollama" else "http://localhost:8000")

    index = []
    queued = []
    tally = collections.Counter()
    for accession, gene in wanted:
        if deadline and time.time() >= deadline:
            log("  time limit reached; run the same command again to continue")
            break
        names = uniprot_names(accession, gene)
        if not names:
            names = [gene]
        classes = classes_for(accession, by_acc)
        query = build_query(cfg, names, classes)
        if args.dry_run:
            log("\n=== %s (%s) classes=%s" % (gene, accession, ",".join(classes) or "-"))
            log("  names: %s" % "; ".join(names))
            log("  query: %s" % query)
            tally["dry"] += 1
            continue

        cache = os.path.join(SEARCH_DIR, safe(gene) + ".json")
        cached = None if args.refresh else read_json(cache)
        if cached:
            hits, results = cached["hitCount"], cached["results"]
        else:
            try:
                hits, results = search(query, pool, cfg)
            except Exception as exc:
                log("  %-10s search FAILED: %s" % (gene, exc))
                tally["search_failed"] += 1
                continue
            write_json(cache, {"gene": gene, "accession": accession, "query": query, "hitCount": hits,
                               "results": results, "retrieved_on": datetime.now().strftime("%Y-%m-%d")})
            time.sleep(PAUSE[0])
        log("\n=== %s: %d hits, considering %d" % (gene, hits, len(results)))
        if not results:
            log("  nothing found. For an antibody-only target (CD19, CD20, BCMA, PSMA, CD3) that is the "
                "correct answer, not a failure: there is no small-molecule assay literature to find.")
            tally["empty"] += 1
            continue

        # Is this protein the subject of the papers, or somebody else's biomarker? Judged by
        # scripts/target_subject.py, which also prunes family-name aliases such as UniProt's
        # "Purinergic receptor" under P2RX3. Without this, the PRMT5 literature queues as MTAP
        # extraction work and files PRMT5 methods under MTAP.
        sub_aliases = tsub.specific_aliases(list(names), gene)
        dropped = [a for a in names if a not in sub_aliases]
        if dropped:
            log("  family alias(es) not used for subject tests: %s" % "; ".join(dropped))
        sub_pat = tsub.alias_pattern(sub_aliases)

        scored = []
        for r in results:
            pmcid = r.get("pmcid")
            if not pmcid:
                continue
            try:
                xml = fetch_fulltext(pmcid)
            except Exception as exc:
                log("  %s full text FAILED: %s" % (pmcid, exc))
                continue
            text, root = body_text(xml)
            if len(text) < 2000:
                continue
            total, hitmap = score_text(text, sig, pen_re, pen_w, small_re)
            verdict = tsub.classify_title(r.get("title") or "", sub_aliases, gene=gene)
            hitmap["title_verdict"] = verdict
            if verdict == "context":
                total += CONTEXT_TITLE_PENALTY
                hitmap["context_title_penalty"] = CONTEXT_TITLE_PENALTY
            scored.append({"pmcid": pmcid, "pmid": r.get("pmid"), "doi": r.get("doi"),
                           "title_verdict": verdict,
                           "title": " ".join((r.get("title") or "").split()),
                           "journal": (r.get("journalInfo") or {}).get("journal", {}).get("title"),
                           "year": r.get("pubYear"), "license": r.get("license"),
                           "score": total, "signals": hitmap, "chars": len(text), "root": root})
        scored.sort(key=lambda x: -x["score"])
        chosen = scored[:keep]
        log("  scored %d papers with full text; keeping %d" % (len(scored), len(chosen)))

        for p in chosen:
            blocks = assay_blocks(p.pop("root")) if p.get("root") is not None else []

            # Which methods blocks name this protein as an assay subject? Judged per block, not per
            # paper: a PRMT5 paper mentions MTAP throughout its introduction, so a whole-document
            # test passes everything, while its methods blocks name MTAP only as the genotype of a
            # cell line. Measured over the MTAP, WRN and TAAR1 caches, this dropped MTAP's four
            # top-scoring papers, all of them PRMT5 work, and kept 20 of 24 WRN papers.
            flags = []
            for name, btext in blocks:
                norm = tsub.normalise(name + " " + btext)
                n_all = len(sub_pat.findall(norm)) if sub_pat else 0
                flags.append("subject" if tsub.subject_positions(norm, sub_pat)
                             else "biomarker_only" if n_all else "no_mention")
            n_subject_blocks = flags.count("subject")
            # Subject blocks first, so --max-blocks truncates the irrelevant ones, not the assay.
            order = {"subject": 0, "no_mention": 1, "biomarker_only": 2}
            paired = sorted(zip(blocks, flags), key=lambda bf: order[bf[1]])
            if args.max_blocks:
                paired = paired[:args.max_blocks]
            blocks = [b for b, _ in paired]
            flags = [f for _, f in paired]
            p["n_blocks"] = len(blocks)
            p["gene"] = gene
            p["subject_blocks"] = n_subject_blocks
            p["biomarker_only_blocks"] = flags.count("biomarker_only")

            # The gate: no methods block names this protein as an assay subject, AND the title does
            # not claim it either. Both conditions, because a block can name the protein once in
            # "Protein expression and purification" and say "the enzyme" thereafter, and a structure
            # paper was dropped by the block test alone during validation.
            if n_subject_blocks == 0 and p["title_verdict"] != "subject":
                log("   [%d] %-9s %-7s  not the subject (%s title, %d/%d blocks name it)  %s"
                    % (p["score"], p["pmcid"], "skip", p["title_verdict"],
                       flags.count("biomarker_only"), len(flags), p["title"][:46]))
                tally["not_subject"] += 1
                continue
            mark = "extract" if p["score"] >= min_score else "skip"
            log("   [%d] %-9s %-7s %2d blocks (%d name the target)  %s"
                % (p["score"], p["pmcid"], mark, len(blocks), n_subject_blocks, p["title"][:52]))
            if p["score"] < min_score:
                tally["below_threshold"] += 1
                continue
            bdir = os.path.join(BLOCK_DIR, p["pmcid"])
            os.makedirs(bdir, exist_ok=True)
            p["blocks"] = []
            for n, (section, text) in enumerate(blocks, 1):
                bpath = os.path.join(bdir, "%02d_%s.txt" % (n, safe(section).replace(" ", "_")[:40]))
                with open(bpath, "w", encoding="utf-8") as handle:
                    handle.write(text)
                p["blocks"].append({"section": section, "path": os.path.relpath(bpath, ROOT),
                                    "chars": len(text), "target": flags[n - 1]})
            tally["papers_kept"] += 1
            tally["blocks"] += len(blocks)

            if args.extract and blocks and args.backend == "claude":
                queued += queue_for_agent(p, gene, blocks)
                tally["queued_blocks"] += len(blocks)
                p.pop("blocks", None)
                continue

            if args.extract and blocks:
                record = {"schema_version": "1.0",
                          "paper": {"pmid": p.get("pmid"), "pmcid": p["pmcid"], "doi": p.get("doi"),
                                    "title": p["title"], "journal": p.get("journal"),
                                    "year": int(p["year"]) if str(p.get("year") or "").isdigit() else None,
                                    "license": p.get("license"), "source": "europepmc_fulltext",
                                    "retrieved_on": datetime.now().strftime("%Y-%m-%d"),
                                    "target_genes": [gene], "paper_type": None, "modality": [],
                                    "relevance_score": min(6, p["score"]),
                                    "extraction": {"model": args.model,
                                                   "run_on": datetime.now().isoformat(timespec="seconds"),
                                                   "passes": ["segment", "extract"], "notes": None}},
                          "assays": [], "unmapped_terms": []}
                for n, (section, text) in enumerate(blocks, 1):
                    if deadline and time.time() >= deadline:
                        log("     time limit reached mid-paper; rerun to continue")
                        break
                    assay_id = "%s-A%d" % (p["pmcid"], n)
                    prompt = build_prompt(p["pmcid"], section, text, assay_id, fn_terms, tech_terms,
                                          args.block_chars)
                    t0 = time.time()
                    try:
                        if args.backend == "ollama":
                            raw = call_ollama(host, args.model, prompt, args.llm_timeout, args.num_ctx)
                        else:
                            raw = call_openai_compatible(host, args.model, prompt, args.llm_timeout, args.api_key)
                    except Exception as exc:
                        log("     %s LLM FAILED: %s" % (assay_id, exc))
                        tally["llm_failed"] += 1
                        continue
                    obj, err = parse_model_json(raw)
                    if obj is None:
                        log("     %s unparseable JSON (%s)" % (assay_id, err))
                        tally["llm_unparseable"] += 1
                        continue
                    obj.setdefault("assay_id", assay_id)
                    obj["assay_id"] = assay_id
                    obj.setdefault("source_section", section)
                    for u in obj.pop("unmapped_terms", []) or []:
                        if isinstance(u, dict):
                            u.setdefault("assay_id", assay_id)
                            record["unmapped_terms"].append(u)
                    record["assays"].append(obj)
                    log("     %s ok in %.0fs (%s / %s)" % (assay_id, time.time() - t0,
                        ",".join(obj.get("technology_labels") or []) or "no tech",
                        obj.get("effect_direction") or "no direction"))
                    tally["assays_extracted"] += 1
                opath = os.path.join(OUT_DIR, "%s.json" % p["pmcid"])
                write_json(opath, record)
                ok, msg = validate(record)
                log("     -> %s (%s)" % (os.path.relpath(opath, ROOT), msg))
                tally["valid" if ok else "invalid"] += 1
            p.pop("blocks", None)
        index += chosen

    if args.dry_run:
        log("\n%d query/queries printed. Nothing fetched." % tally["dry"])
        return
    for p in index:
        p.pop("root", None)
    write_json(os.path.join(OUT_DIR, "index.json"),
               {"built": datetime.now().isoformat(timespec="seconds"), "papers": index})
    log("\n--- summary ---")
    for k in ("papers_kept", "blocks", "not_subject", "below_threshold", "empty", "search_failed",
              "assays_extracted", "queued_blocks", "valid", "invalid", "llm_failed", "llm_unparseable"):
        if tally[k]:
            log("  %-18s %d" % (k, tally[k]))
    log("  index: %s" % os.path.relpath(os.path.join(OUT_DIR, "index.json"), ROOT))

    if queued:
        tasks = write_agent_tasks(queued)
        todo = sum(1 for t in queued if not t["done"])
        log("\nQueued %d assay block(s) for the Claude Code agent (%d still to do)." % (len(queued), todo))
        log("  %s" % os.path.relpath(tasks, ROOT))
        log("\nIn Claude Code, say:")
        log('  "Work the queue in build_output/papers/agent_queue/TASKS.md"')
        log("Then merge and validate what the agent wrote:")
        log("  py scripts\\fetch_papers.py --agent-merge")
    elif not args.extract and tally["papers_kept"]:
        log("\nBlocks are in build_cache/papers/blocks/. To fill the JSON, either:")
        log("  # a local model (cheap, needs tuning)")
        log("  py scripts\\fetch_papers.py %s --extract --model qwen2.5:7b-instruct"
            % ("--genes " + args.genes if args.genes else "--all"))
        log("  # or the Claude Code agent in this session (better, for a handful of papers)")
        log("  py scripts\\fetch_papers.py %s --extract --backend claude"
            % ("--genes " + args.genes if args.genes else "--all"))


# ------------------------------------------------------------------ the agent path (--backend claude)
#
# The script cannot call Claude Code: Claude Code is the agent running the script. So this backend inverts
# the direction. The script does the parts it is good at (search, fetch, score, split, validate, merge) and
# leaves a work queue; the agent reads each block and writes one JSON file per assay; --agent-merge then
# assembles and validates. Use it for a handful of papers where the extraction quality matters more than
# the throughput - a frontier model reading 45 nested fields out of a methods paragraph will beat a 7B,
# but it is the expensive path, so it is not the one to point at 5,852 proteins.

AGENT_DIR = os.path.join(OUT_DIR, "agent_queue")


def queue_for_agent(paper, gene, blocks):
    """Write the queue entry for one paper: the paper metadata the agent should not have to retype,
    and one todo per assay block. Returns the todo list."""
    pmcid = paper["pmcid"]
    pdir = os.path.join(AGENT_DIR, pmcid)
    os.makedirs(pdir, exist_ok=True)
    todos = []
    for n, (section, _text) in enumerate(blocks, 1):
        assay_id = "%s-A%d" % (pmcid, n)
        todos.append({
            "assay_id": assay_id,
            "pmcid": pmcid,
            "gene": gene,
            "section": section,
            "block_path": paper["blocks"][n - 1]["path"],
            "block_chars": paper["blocks"][n - 1]["chars"],
            "write_to": os.path.relpath(os.path.join(pdir, "A%d.json" % n), ROOT),
            "done": os.path.exists(os.path.join(pdir, "A%d.json" % n)),
        })
    write_json(os.path.join(pdir, "paper.json"), {
        "pmid": paper.get("pmid"), "pmcid": pmcid, "doi": paper.get("doi"),
        "title": paper["title"], "journal": paper.get("journal"),
        "year": int(paper["year"]) if str(paper.get("year") or "").isdigit() else None,
        "license": paper.get("license"), "source": "europepmc_fulltext",
        "retrieved_on": datetime.now().strftime("%Y-%m-%d"),
        "target_genes": [gene], "relevance_score": min(6, paper["score"]),
        "score_signals": paper.get("signals") or {},
    })
    return todos


AGENT_TASKS_HEADER = """# Assay extraction queue for the Claude Code agent

Generated by `scripts/fetch_papers.py --extract --backend claude` on __WHEN__.
__NTODO__ assay block(s) across __NPAPER__ paper(s). __NDONE__ already done.

## What to do

For each row in the table below that is not done:

1. Read the block file with the Read tool. It is one methods subsection: one assay.
2. Write ONE json file to `write_to`, holding ONE assay object that conforms to the `assay`
   definition in `schema/paper_extraction.schema.json`. Do not wrap it in a paper record and do not
   put several assays in one file - `--agent-merge` assembles those.
3. When every block of a paper is written, run:

   ```
   py scripts\\fetch_papers.py --agent-merge
   ```

   That builds `build_output/papers/<pmcid>.json` from `paper.json` plus the A*.json files, hoists
   `unmapped_terms` to the record, and validates the whole thing against the schema. It reports any
   error with its field path. Fix and rerun; it is idempotent.

A block already written is skipped, so this is resumable: do a few, merge, come back.

## The rules that matter

These are the same rules the local-model prompt carries, and they are where extractions go wrong:

1. **Never guess.** A field the block does not state is `null`, AND its dotted path goes in
   `not_reported`. A null with no `not_reported` entry means "missed"; the difference is the whole
   point of the field.
2. **Final, not stock.** "compounds dissolved in 100% DMSO" is a stock, so `conditions.dmso_pct`
   stays null. "1% DMSO final" gives `1`.
3. **Do not borrow from a neighbouring assay.** Papers describe several assays and often share a
   buffer paragraph. Only record what belongs to THIS block. If the block relies on a buffer defined
   in a different section, record it and add `method_by_reference` to `extraction_flags`.
4. **Residue ranges need a cue.** "residues 500-942" or "helicase core (500-942)" is a construct.
   A bare "(2015-2018)" or "(12-17)" is a citation year or page range: leave `construct_start` null.
5. **Activators get the same care as inhibitors.** If the assay reads activation, set
   `effect_direction` to `activation` or `both`, fill `pharmacology`, and fill `stimulation` with the
   reference agonist and its level (EC20 / EC50 / EC80 / ECmax / unstimulated). An EC50 from a
   potentiator screen is meaningless without the challenge level. Fill `protocol_additions` for a
   multi-add protocol, and `controls.basal` / `controls.max_stimulation` rather than forcing an
   activation assay into positive/negative controls.
6. **Quote exactly.** `evidence` needs at least one verbatim quote copied from the block, with the
   field path it supports. This is what makes the record reviewable without reopening the paper.
7. **Do not invent vocabulary.** `function_labels` and `technology_labels` must come from
   `scripts/term_table.csv` (55 function terms, 32 technology terms). Anything you cannot map goes in
   this assay's `unmapped_terms` array with its quote, not into the nearest enum.
8. A DMSO curve run to correct an instrument artefact (SPR refractive-index matching) is NOT a DMSO
   tolerance test. Leave `quality.dmso_tolerance_tested_pct` empty and say so in the notes.

`schema/paper_extraction.example.json` is a filled inhibitor paper and
`schema/paper_extraction.example-activator.json` a filled activator paper. Copy their shape.

## Queue

| done | assay_id | gene | chars | section | block_path | write_to |
|---|---|---|---|---|---|---|
"""


def write_agent_tasks(todos):
    os.makedirs(AGENT_DIR, exist_ok=True)
    papers = sorted({t["pmcid"] for t in todos})
    ndone = sum(1 for t in todos if t["done"])
    body = AGENT_TASKS_HEADER
    for token, value in (("__WHEN__", datetime.now().strftime("%Y-%m-%d %H:%M")),
                         ("__NTODO__", str(len(todos))),
                         ("__NPAPER__", str(len(papers))),
                         ("__NDONE__", str(ndone))):
        body = body.replace(token, value)
    rows = []
    for t in todos:
        rows.append("| %s | `%s` | %s | %d | %s | `%s` | `%s` |" % (
            "x" if t["done"] else " ", t["assay_id"], t["gene"], t["block_chars"],
            t["section"].replace("|", "/")[:54], t["block_path"], t["write_to"]))
    path = os.path.join(AGENT_DIR, "TASKS.md")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(body + "\n".join(rows) + "\n")
    write_json(os.path.join(AGENT_DIR, "queue.json"),
               {"built": datetime.now().isoformat(timespec="seconds"), "todos": todos})
    return path


def agent_merge():
    """Assemble the agent's per-assay JSON files into one record per paper, then validate."""
    if not os.path.isdir(AGENT_DIR):
        log("No agent queue at %s. Run with --extract --backend claude first."
            % os.path.relpath(AGENT_DIR, ROOT))
        return 1
    n_ok = n_bad = n_waiting = 0
    for pmcid in sorted(os.listdir(AGENT_DIR)):
        pdir = os.path.join(AGENT_DIR, pmcid)
        if not os.path.isdir(pdir):
            continue
        meta = read_json(os.path.join(pdir, "paper.json"))
        if not meta:
            log("  %-12s no paper.json, skipped" % pmcid)
            continue
        files = sorted(f for f in os.listdir(pdir)
                       if re.fullmatch(r"A\d+\.json", f))
        if not files:
            n_waiting += 1
            log("  %-12s no assay files written yet" % pmcid)
            continue
        assays, unmapped, bad = [], [], []
        for f in files:
            obj = read_json(os.path.join(pdir, f))
            if not isinstance(obj, dict):
                bad.append("%s is not a JSON object" % f)
                continue
            if "assays" in obj or "paper" in obj:
                bad.append("%s looks like a whole paper record; it must hold ONE assay object" % f)
                continue
            aid = "%s-A%s" % (pmcid, re.sub(r"\D", "", f))
            obj.setdefault("assay_id", aid)
            for u in obj.pop("unmapped_terms", []) or []:
                if isinstance(u, dict):
                    u.setdefault("assay_id", obj["assay_id"])
                    unmapped.append(u)
            assays.append(obj)
        record = {"schema_version": "1.0",
                  "paper": {k: meta.get(k) for k in ("pmid", "pmcid", "doi", "title", "journal", "year",
                                                     "license", "source", "retrieved_on", "target_genes",
                                                     "relevance_score")},
                  "assays": assays, "unmapped_terms": unmapped}
        record["paper"]["paper_type"] = meta.get("paper_type")
        record["paper"]["modality"] = meta.get("modality") or []
        record["paper"]["extraction"] = {"model": "claude-code-agent",
                                         "run_on": datetime.now().isoformat(timespec="seconds"),
                                         "passes": ["segment", "extract-by-agent", "merge"],
                                         "notes": meta.get("notes")}
        opath = os.path.join(OUT_DIR, "%s.json" % pmcid)
        write_json(opath, record)
        ok, msg = validate(record)
        for b in bad:
            log("  %-12s WARNING %s" % (pmcid, b))
        log("  %-12s %d assay(s) -> %s (%s)" % (pmcid, len(assays), os.path.relpath(opath, ROOT), msg))
        n_ok += 1 if ok and not bad else 0
        n_bad += 0 if ok and not bad else 1
    log("\nmerged: %d clean, %d with problems, %d still waiting for the agent" % (n_ok, n_bad, n_waiting))
    return 0 if n_bad == 0 else 2


def validate_path(path):
    """Validate one record, or every *.json in a directory."""
    target = path if os.path.isabs(path) else os.path.join(ROOT, path)
    files = []
    if os.path.isdir(target):
        files = [os.path.join(target, f) for f in sorted(os.listdir(target))
                 if f.endswith(".json") and f != "index.json"]
    elif os.path.exists(target):
        files = [target]
    else:
        log("Not found: %s" % target)
        return 1
    bad = 0
    for f in files:
        rec = read_json(f)
        if rec is None:
            log("  %-34s UNREADABLE / not JSON" % os.path.basename(f))
            bad += 1
            continue
        ok, msg = validate(rec)
        log("  %-34s %s" % (os.path.basename(f), msg))
        bad += 0 if ok else 1
    log("\n%d file(s), %d with schema errors" % (len(files), bad))
    return 0 if bad == 0 else 2


def validate(record):
    try:
        import jsonschema
    except ImportError:
        return True, "not validated (pip install jsonschema)"
    schema = read_json(SCHEMA)
    if not schema:
        return True, "not validated (schema unreadable)"
    errs = sorted(jsonschema.Draft202012Validator(schema).iter_errors(record),
                  key=lambda e: list(e.path))
    if not errs:
        return True, "%d assay(s), valid against the schema" % len(record["assays"])
    first = errs[0]
    return False, "%d schema error(s); first at /%s: %s" % (
        len(errs), "/".join(str(x) for x in first.path), first.message[:120])


if __name__ == "__main__":
    main()
