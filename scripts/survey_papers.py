#!/usr/bin/env python
"""Survey how much open-access assay literature exists per protein, and rank the proteins by how much
a paper-mining pass would actually add to the database.

Two halves:

  demand   Local only, no network. For each protein the app knows, how thin is its CURRENT data?
           A protein with 400 campaigns and full condition coverage has little to gain from papers;
           a real drug target sitting below the app's own thresholds has a lot. The thresholds are
           the app's, not invented here: MINC = 5 classifiable campaigns before the website shows
           gene-level usage rather than falling back to the function class, and COND_MIN_N = 10
           campaigns stating a parameter before a gene-level condition row is written.

  supply   One lightweight Europe PMC call per protein: hitCount plus the first page of titles.
           No full text is fetched. Uses the same query shape as fetch_papers.py (target in
           TITLE_ABS, assay terms in METHODS), so the count is of papers whose methods section
           really does describe an assay. Discounted by scripts/target_subject.py, which asks
           whether the protein is the subject of those papers or just the biomarker that made
           somebody else's target interesting.

The score multiplies the two, because either alone is misleading: papers that do not exist cannot be
mined, and papers about an already-covered protein add nothing.

Examples
--------
  # demand side only, instant, writes the worklist with supply unknown
  py scripts\\survey_papers.py --demand-only

  # survey the proteins with at least 10 campaigns (about 1,500 calls)
  py scripts\\survey_papers.py --min-campaigns 10

  # everything with any campaign at all, resumable; run it in the background
  py scripts\\survey_papers.py --min-campaigns 1 --max-minutes 120

  # rebuild the ranking from the cache without any new calls
  py scripts\\survey_papers.py --rank-only

Outputs
-------
  build_output/papers/survey.csv    one row per protein, every input to the score, sortable
  build_output/papers/survey.md     the prioritised shortlist, grouped, with the reasoning
  build_cache/papers/survey/<GENE>.json   cached hitCount + titles, so reruns cost nothing
"""

import argparse
import collections
import csv
import glob
import json
import math
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
CACHE = os.path.join(ROOT, "build_cache")
SURVEY_DIR = os.path.join(CACHE, "papers", "survey")
OUT_DIR = os.path.join(ROOT, "build_output", "papers")
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA = "bioassay-builder literature survey (research use)"

MINC = 5          # docs/index.html: fewer classifiable campaigns than this and the site uses the class
COND_MIN_N = 10   # scripts/build_database.py COND_MIN_N["gene"]

sys.path.insert(0, HERE)
import target_subject as tsub
try:
    import fetch_papers as fp
except ImportError:
    fp = None


def log(msg):
    print(msg)
    sys.stdout.flush()


def read_json(path, default=None):
    try:
        with open(path, encoding="utf-8") as h:
            return json.load(h)
    except (OSError, ValueError):
        return default


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as h:
        json.dump(obj, h, indent=1, ensure_ascii=False)


def read_csv_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))


# ------------------------------------------------------------------ demand: what the app lacks today

def build_demand():
    """One record per protein with a cached summary. Local files only."""
    assays = read_csv_rows(os.path.join(DATA, "assays.csv"))
    entry_classes = {}
    entry_labels = {}
    for a in assays:
        eid = a["entry_id"]
        entry_classes[eid] = {c.strip() for c in (a.get("classes") or "").split(";") if c.strip()}
        labs = set()
        for part in (a.get("evidence_labels") or a.get("evidence_function") or "").split(";"):
            for x in part.split("|"):
                if x.strip():
                    labs.add(x.strip())
        entry_labels[eid] = labs

    # gene-level condition coverage already achieved
    cond_params = collections.defaultdict(set)
    for r in read_csv_rows(os.path.join(DATA, "condition_stats.csv")):
        if (r.get("scope_type") or "") == "gene" and r.get("scope_id"):
            cond_params[r["scope_id"]].add(r.get("parameter") or "")

    # gene-level usage rows already achieved, and which readouts are known per gene
    usage_genes = collections.defaultdict(int)
    usage_readouts = collections.defaultdict(set)
    for r in read_csv_rows(os.path.join(DATA, "usage_stats.csv")):
        if (r.get("scope_type") or "") != "gene" or not r.get("scope_id"):
            continue
        g = r["scope_id"]
        rd = r.get("readout") or ""
        if not rd:
            usage_genes[g] += 1
        elif not rd.startswith("plate:"):
            usage_readouts[g].add(rd)

    rows = []
    for path in sorted(glob.glob(os.path.join(CACHE, "summary", "*.json"))):
        s = read_json(path)
        if not s:
            continue
        gene = (s.get("scope") or s.get("gene") or "").strip()
        if not gene:
            continue
        camps = s.get("campaigns") or []
        classifiable = [c for c in camps if c.get("f")]
        funcs = collections.Counter()
        techs = collections.Counter()
        n_state_tech = 0
        for c in classifiable:
            for f in c.get("f") or []:
                funcs[f] += 1
            t = c.get("t") or []
            if t:
                n_state_tech += 1
            for x in t:
                techs[x] += 1
        classes = s.get("classes") or []
        # how many playbook entries could ever apply to this protein
        n_entries = sum(1 for eid, cl in entry_classes.items() if cl & (set(classes) | {"generic"}))
        # of those, how many have an evidence label this protein's campaigns actually measured
        measured = {f for f in funcs}
        n_entries_with_data = sum(1 for eid, cl in entry_classes.items()
                                  if (cl & (set(classes) | {"generic"})) and (entry_labels[eid] & measured))
        rows.append({
            "gene": gene,
            "accession": s.get("acc") or "",
            "classes": ";".join(classes),
            "n_classes": len(classes),
            "n_campaigns": len(camps),
            "n_classifiable": len(classifiable),
            "n_stating_technology": n_state_tech,
            "n_distinct_readouts": len(techs),
            "n_condition_params": len(cond_params.get(gene, set())),
            "n_usage_rows": usage_genes.get(gene, 0),
            "n_entries_possible": n_entries,
            "n_entries_with_data": n_entries_with_data,
            "below_minc": 1 if len(classifiable) < MINC else 0,
            "below_cond_min": 1 if len(cond_params.get(gene, set())) == 0 else 0,
            "top_function": (funcs.most_common(1) or [("", 0)])[0][0],
        })
    return rows


def demand_score(r):
    """How much this protein stands to gain. Higher = thinner current data on a protein that is
    nonetheless really worked on."""
    n = r["n_classifiable"]
    # Real target: campaigns exist at all. Flat above ~50, because a protein with 400 campaigns is
    # not 8x more worth mining than one with 50 - the app already has plenty for it.
    worked_on = min(1.0, math.log10(n + 1) / math.log10(51))
    # Gaps, each 0..1
    gap_usage = 1.0 if n < MINC else 0.0                       # site falls back to the class
    gap_cond = 1.0 if r["n_condition_params"] == 0 else max(0.0, 1 - r["n_condition_params"] / 8.0)
    gap_readout = 1.0 if r["n_distinct_readouts"] == 0 else max(0.0, 1 - r["n_distinct_readouts"] / 6.0)
    gap_entries = 0.0
    if r["n_entries_possible"]:
        gap_entries = 1 - r["n_entries_with_data"] / r["n_entries_possible"]
    gap = 0.35 * gap_cond + 0.25 * gap_readout + 0.25 * gap_entries + 0.15 * gap_usage
    return round(worked_on * gap, 4), round(worked_on, 4), round(gap, 4)


# ------------------------------------------------------------------ supply: one light EPMC call each

def survey_query(cfg, names, class_ids):
    terms = list(cfg["assay_terms_core"])
    for c in class_ids:
        terms += cfg["assay_terms_by_class"].get(c, [])
    seen, uniq = set(), []
    for t in terms:
        if t not in seen:
            seen.add(t)
            uniq.append(t)
    name_q = " OR ".join('TITLE_ABS:"%s"' % n for n in names)
    term_q = " OR ".join("METHODS:%s" % t for t in uniq)
    return cfg["query_template"].format(target_names=name_q, assay_terms=term_q,
                                        exclusions=cfg.get("exclusions", "")).strip()


def epmc_light(query, n_titles, timeout=60, tries=3):
    """hitCount plus the first page of titles. resultType=lite keeps the payload small; no full text."""
    params = {"query": query, "format": "json", "resultType": "lite", "pageSize": max(1, min(n_titles, 25))}
    url = EPMC + "?" + urllib.parse.urlencode(params)
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                d = json.load(resp)
            res = (d.get("resultList") or {}).get("result") or []
            return {"hitCount": d.get("hitCount", 0),
                    "titles": [{"pmid": x.get("pmid"), "pmcid": x.get("pmcid"),
                                "year": x.get("pubYear"),
                                "title": " ".join((x.get("title") or "").split())[:240]}
                               for x in res]}
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code == 400:
                raise
            time.sleep(2 ** i)
    raise last


# ------------------------------------------------------------------ scoring and output

# Titles that say "this paper develops or characterises an assay" rather than "this paper is about biology"
ASSAY_TITLE_RE = re.compile(
    r"\bassay\b|high[- ]throughput|\bhts\b|\bscreen(ing|ed|s)?\b|inhibitor|antagonist|agonist|"
    r"allosteric|potentiator|activator|\bprobe\b|fragment|\bsar\b|discovery of|identification of|"
    r"optimi[sz]ation|structure[- ]activity|biochemical|biophysical|binding|potency|\bic50\b|"
    r"\bec50\b|crystal structure|co-?crystal|degrader|\bprotac\b", re.I)
# A genetic screen is not an assay-development paper, but it matches "screen". Checked first, and it
# cancels the assay match: "CRISPR Screening Identifies ... Modulators of WRN Dependency" is target
# biology, however much assay vocabulary the title carries.
GENETIC_SCREEN_RE = re.compile(
    r"crispr|\brnai\b|\bsirna\b|shrna|genome[- ]wide|genetic (screen|dependenc)|"
    r"knockout screen|dependency map|\bdepmap\b", re.I)
BIOLOGY_TITLE_RE = re.compile(
    r"\bin (male |female |adult |young |aged )?mice\b|\bmouse model\b|\bknockout\b|\bcohort\b|"
    r"\bpatients?\b|prognos|clinical outcome|meta-analysis|expression (in|profile|level)|"
    r"single[- ]cell|transcriptom|epidemiolog|case report|survival|biomarker of|immunohistochem|"
    r"xenograft|\bin vivo\b|pharmacokinetic|neuroinflamm|attenuates|ameliorates|alleviates|"
    r"suppresses .* (progression|growth|tumor)|promotes|induces|modulat\w+ .* (axis|pathway|signaling)",
    re.I)


def title_quality(titles, verdicts=None):
    """Of the titles we got back, what fraction look like assay / discovery work rather than biology?
    A cheap stand-in for the full-text score, using only what a lightweight query returns.

    `verdicts` is the parallel output of target_subject.classify_title. A title whose only mention of
    the protein is a biomarker mention can never count as assay work for THIS protein, however much
    assay vocabulary it carries: "Discovery of selective PRMT5 inhibitors for MTAP-deleted cancers"
    is assay work for PRMT5, and counting it for MTAP is what put MTAP fourth out of 3,441."""
    if not titles:
        return 0.0, 0, 0
    assayish = biologyish = 0
    for n, t in enumerate(titles):
        s = t.get("title") or ""
        if verdicts and verdicts[n] == "context":
            biologyish += 1
            continue
        a = bool(ASSAY_TITLE_RE.search(s)) and not GENETIC_SCREEN_RE.search(s)
        b = bool(BIOLOGY_TITLE_RE.search(s)) or bool(GENETIC_SCREEN_RE.search(s))
        if a and not b:
            assayish += 1
        elif b and not a:
            biologyish += 1
    # Titles that match both, or neither, count in the denominator but not the numerator, so an
    # ambiguous page dilutes the score rather than inflating it.
    return round(assayish / len(titles), 3), assayish, biologyish


def supply_score(hits, frac_assay, frac_subject=1.0):
    """Papers available, damped: 0 papers is useless, 30 is plenty, 3,000 is not 100x better than 30.

    frac_subject is the final guard. Without it the 0.35 floor below gives a protein credit for a
    large literature even when none of that literature assays the protein; with it, a gene that only
    ever appears as somebody else's stratification marker cannot rank on volume alone."""
    if not hits:
        return 0.0
    volume = min(1.0, math.log10(hits + 1) / math.log10(31))
    # a big pile of biology papers is worth less than a small pile of assay papers
    return round(volume * (0.35 + 0.65 * frac_assay) * frac_subject, 4)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-campaigns", type=int, default=10,
                    help="only survey proteins with at least this many classifiable campaigns (default 10, "
                         "about 1,500 proteins; 1 covers about 5,140)")
    ap.add_argument("--limit", type=int, default=0, help="stop after this many API calls")
    ap.add_argument("--titles", type=int, default=10, help="titles to keep per protein (default 10)")
    ap.add_argument("--pause", type=float, default=1.0, help="seconds between calls")
    ap.add_argument("--max-minutes", type=float, default=0, help="stop after this long (resumable)")
    ap.add_argument("--refresh", action="store_true", help="re-query proteins already cached")
    ap.add_argument("--demand-only", action="store_true", help="no network; write the demand side only")
    ap.add_argument("--rank-only", action="store_true", help="no network; rebuild outputs from the cache")
    ap.add_argument("--top", type=int, default=120, help="rows in the markdown shortlist")
    args = ap.parse_args(argv)

    os.makedirs(SURVEY_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    log("Building the demand side from local files...")
    rows = build_demand()
    for r in rows:
        r["demand"], r["worked_on"], r["gap"] = demand_score(r)
    log("  %d proteins; %d with >=1 classifiable campaign; %d with >=%d"
        % (len(rows), sum(1 for r in rows if r["n_classifiable"] >= 1),
           sum(1 for r in rows if r["n_classifiable"] >= args.min_campaigns), args.min_campaigns))

    by_gene = {r["gene"]: r for r in rows}
    cfg = None
    if not args.demand_only:
        cfg = read_json(os.path.join(DATA, "paper_queries.json"))
        if not cfg:
            sys.exit("Cannot read data/paper_queries.json")

    # Query in demand order, so an interrupted run has done the proteins that matter most.
    todo = [r for r in rows if r["n_classifiable"] >= args.min_campaigns]
    todo.sort(key=lambda r: -r["demand"])
    n_call = n_cached = 0
    deadline = time.time() + args.max_minutes * 60 if args.max_minutes else None

    if not args.demand_only:
        log("\nSurveying Europe PMC (hitCount + %d titles, no full text)..." % args.titles)
        for n, r in enumerate(todo, 1):
            if args.limit and n_call >= args.limit:
                log("  call limit reached")
                break
            if deadline and time.time() >= deadline:
                log("  time limit reached; rerun the same command to continue")
                break
            cpath = os.path.join(SURVEY_DIR, re.sub(r"[^A-Za-z0-9_.-]", "_", r["gene"]) + ".json")
            cached = None if args.refresh else read_json(cpath)
            if cached and not args.rank_only:
                n_cached += 1
                continue
            if args.rank_only:
                continue
            names = fp.uniprot_names(r["accession"], r["gene"]) if fp else [r["gene"]]
            if not names:
                names = [r["gene"]]
            classes = [c for c in (r["classes"] or "").split(";") if c]
            q = survey_query(cfg, names, classes)
            try:
                res = epmc_light(q, args.titles)
            except Exception as exc:
                log("  %-12s FAILED %s" % (r["gene"], exc))
                write_json(cpath, {"gene": r["gene"], "error": str(exc),
                                   "surveyed_on": datetime.now().strftime("%Y-%m-%d")})
                n_call += 1
                time.sleep(args.pause)
                continue
            res.update({"gene": r["gene"], "accession": r["accession"], "query": q,
                        "names": names, "surveyed_on": datetime.now().strftime("%Y-%m-%d")})
            write_json(cpath, res)
            n_call += 1
            if n_call % 50 == 0:
                log("  %d called, %d cached, at %s (%d hits)" % (n_call, n_cached, r["gene"], res["hitCount"]))
            time.sleep(args.pause)
        log("  %d new call(s), %d already cached" % (n_call, n_cached))

    # fold the cache in
    n_have = 0
    for path in glob.glob(os.path.join(SURVEY_DIR, "*.json")):
        d = read_json(path) or {}
        g = d.get("gene")
        if g not in by_gene:
            continue
        r = by_gene[g]
        if d.get("error"):
            r["epmc_hits"] = ""
            r["survey_error"] = d["error"][:60]
            continue
        titles = d.get("titles") or []
        # Is this protein the subject of its own literature, or somebody else's biomarker?
        names = d.get("names") or [g]
        aliases = tsub.specific_aliases(list(names), g)
        apat = tsub.alias_pattern(aliases)
        verdicts = [tsub.classify_title(t.get("title"), aliases, apat, g) for t in titles]
        subj_frac, n_subj, n_ctx, n_abs = tsub.subject_fraction(titles, names, gene=g)
        frac, assayish, biologyish = title_quality(titles, verdicts)
        r["epmc_hits"] = d.get("hitCount", 0)
        r["titles_seen"] = len(titles)
        r["titles_assayish"] = assayish
        r["titles_biologyish"] = biologyish
        r["frac_assay_titles"] = frac
        r["titles_subject"] = n_subj
        r["titles_context"] = n_ctx
        r["titles_no_name"] = n_abs
        r["frac_subject"] = subj_frac
        r["aliases_dropped"] = ";".join(a for a in names if a not in aliases)
        r["supply"] = supply_score(d.get("hitCount", 0), frac, subj_frac)
        r["example_title"] = (titles[0]["title"] if titles else "")[:150]
        n_have += 1
    log("\n%d protein(s) have a survey result" % n_have)

    for r in rows:
        r.setdefault("epmc_hits", "")
        for k in ("titles_seen", "titles_assayish", "titles_biologyish",
                  "titles_subject", "titles_context", "titles_no_name"):
            r.setdefault(k, "")
        r.setdefault("frac_assay_titles", "")
        r.setdefault("frac_subject", "")
        r.setdefault("aliases_dropped", "")
        r.setdefault("supply", "")
        r.setdefault("example_title", "")
        r.setdefault("survey_error", "")
        r["priority"] = round(r["demand"] * r["supply"], 4) if r["supply"] != "" else ""

    cols = ["gene", "accession", "classes", "n_classifiable", "n_stating_technology", "n_distinct_readouts",
            "n_condition_params", "n_usage_rows", "n_entries_possible", "n_entries_with_data",
            "below_minc", "below_cond_min", "top_function", "worked_on", "gap", "demand",
            "epmc_hits", "titles_seen", "titles_assayish", "titles_biologyish", "frac_assay_titles",
            "titles_subject", "titles_context", "titles_no_name", "frac_subject",
            "supply", "priority", "example_title", "aliases_dropped", "survey_error"]
    out_csv = os.path.join(OUT_DIR, "survey.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (-(r["priority"] if r["priority"] != "" else -1), -r["demand"])):
            w.writerow(r)
    log("  %s" % os.path.relpath(out_csv, ROOT))
    write_markdown(rows, args, os.path.join(OUT_DIR, "survey.md"))
    log("  %s" % os.path.relpath(os.path.join(OUT_DIR, "survey.md"), ROOT))
    return 0


def write_markdown(rows, args, path):
    scored = [r for r in rows if r["priority"] != ""]
    scored.sort(key=lambda r: -r["priority"])
    done = len(scored)
    lines = []
    lines.append("# Which proteins to mine papers for, in order\n")
    lines.append("Generated by `scripts/survey_papers.py` on %s. %d protein(s) surveyed "
                 "(Europe PMC hit counts and titles only; no full text fetched).\n"
                 % (datetime.now().strftime("%Y-%m-%d %H:%M"), done))
    lines.append("""## How the order was decided

`priority = demand x supply`, because either alone misleads.

**demand** = `worked_on x gap`, from our own data only.
- `worked_on` = log-damped classifiable campaign count, flat above ~50. A protein nobody screens is
  not worth mining however many papers exist; a protein with 400 campaigns is not eight times more
  worth mining than one with 50, because we already have plenty for it.
- `gap` = what we are missing, weighted: 35 per cent no gene-level condition rows (a protein needs
  __CONDMIN__ campaigns stating a parameter before `condition_stats.csv` gets a row), 25 per cent few
  distinct readouts, 25 per cent playbook entries that apply to this protein but have no evidence for
  it, 15 per cent sitting below MINC = __MINC__ so the website shows class-level numbers instead of
  the protein's own.

**supply** = log-damped Europe PMC hit count, multiplied by the fraction of returned titles that look
like assay or discovery work rather than biology, then by `frac_subject`. Titles are a cheap stand-in
for the full-text score in `data/paper_queries.json`; they are not as good, which is the main caveat
on this ranking.

**frac_subject** (`scripts/target_subject.py`) asks whether the protein is the subject of its own
literature or somebody else's biomarker. "Discovery of selective PRMT5 inhibitors for MTAP-deleted
cancers" is assay work for PRMT5; mining it for MTAP readouts files PRMT5 methods under MTAP. Each
title is judged subject, context or absent from the wording around the protein's name, a context
title can never count as assay work, and the fraction is shrunk toward 1.0 so one biomarker title
does not sink a protein. The same module drops family-name aliases: UniProt lists "Purinergic
receptor" under P2RX3 and "CK1" under CSNK1A1L, and those match every sibling in the family.
`aliases_dropped` in the CSV records what was removed for each protein.
""".replace("__CONDMIN__", str(COND_MIN_N)).replace("__MINC__", str(MINC)))

    def table(rs, n):
        out = ["| # | gene | classes | campaigns | cond params | readouts | entries w/ data | EPMC hits | assay titles | subj/ctx | priority |",
               "|---|---|---|---|---|---|---|---|---|---|---|"]
        for i, r in enumerate(rs[:n], 1):
            out.append("| %d | **%s** | %s | %d | %d | %d | %d/%d | %s | %s/%s | %s/%s | %.3f |" % (
                i, r["gene"], (r["classes"] or "-").replace(";", ", ")[:34], r["n_classifiable"],
                r["n_condition_params"], r["n_distinct_readouts"], r["n_entries_with_data"],
                r["n_entries_possible"], r["epmc_hits"], r["titles_assayish"], r["titles_seen"],
                r["titles_subject"], r["titles_context"], r["priority"]))
        return "\n".join(out)

    lines.append("\n## Top %d\n" % min(args.top, done))
    lines.append(table(scored, args.top))

    # Grouped by what a paper would actually add, which is more actionable than one flat list.
    lines.append("\n\n## The same list, grouped by what mining would add\n")
    groups = [
        ("No conditions of its own, and papers exist",
         "These show class-level buffer / pH / DMSO on the site today. A methods section gives them their own.",
         lambda r: r["n_condition_params"] == 0 and r["n_classifiable"] >= MINC),
        ("Readout list is thin",
         "Few distinct readouts seen in campaigns, so the 'Most used for <protein>' line has little to say.",
         lambda r: r["n_distinct_readouts"] <= 3 and r["n_condition_params"] > 0),
        ("Playbook entries apply but have no evidence for this protein",
         "Entries are shown from the function class with nothing protein-specific behind them.",
         lambda r: r["n_entries_possible"] and r["n_entries_with_data"] / r["n_entries_possible"] < 0.5),
        ("Below MINC, so the site falls back to the class",
         "Fewer than %d classifiable campaigns. Papers are the only route to protein-specific numbers, "
         "but check it is a real programme first." % MINC,
         lambda r: r["n_classifiable"] < MINC),
    ]
    used = set()
    for title, why, test in groups:
        sel = [r for r in scored if test(r) and r["gene"] not in used][:25]
        for r in sel:
            used.add(r["gene"])
        lines.append("\n### %s (%d shown)\n" % (title, len(sel)))
        lines.append("%s\n" % why)
        if sel:
            lines.append(table(sel, 25))
        else:
            lines.append("_none_")

    lines.append("\n\n## Read this before working down the list\n")
    big = [r for r in scored if r["n_classifiable"] >= 50]
    thin = [r for r in scored if r["n_classifiable"] < MINC]
    nocond = [r for r in scored if r["n_condition_params"] == 0]
    nopapers = [r for r in scored if r["epmc_hits"] == 0]
    lines.append("- **%d of the %d surveyed have no gene-level condition rows at all.** That is the single "
                 "biggest gap, and it is why `gap` weights it most: those proteins currently show "
                 "class-level conditions on the site.\n" % (len(nocond), done))
    lines.append("- **%d have zero Europe PMC hits** with the assay-scoped query. Mining cannot help them; "
                 "they are at the bottom of the list and should stay there.\n" % len(nopapers))
    lines.append("- **%d are heavily screened already (>=50 campaigns).** They score lower on `worked_on` "
                 "gain than their fame suggests, deliberately: we already have their usage and "
                 "conditions, so a paper adds construct detail and Z'/DMSO, not usage.\n" % len(big))
    lines.append("- **%d sit below MINC = %d.** The website shows them class-level numbers. Papers are the "
                 "only way to give them anything of their own, but check they are real targets first: "
                 "many are proteins with one or two ChEMBL assays and no programme behind them.\n"
                 % (len(thin), MINC))
    lines.append("\nFull table with every input column: `build_output/papers/survey.csv`.\n")
    with open(path, "w", encoding="utf-8") as h:
        h.write("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
