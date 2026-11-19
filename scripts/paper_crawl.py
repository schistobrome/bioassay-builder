"""Paper crawl: pick targets in batches, read their open-access papers, and log two things -

  1. evidence that an assay type is used for a target   -> crawl/evidence.csv
  2. technologies and assay types our vocabulary lacks  -> crawl/term_occurrences.csv, catalogued
                                                            by hand in crawl/term_catalog.csv

The process, and the prompt that runs it, are in docs/PAPER_CRAWL.md. This script does the
bookkeeping; the agent (Claude Code) does the reading.

Subcommands
-----------
  pick      --n 20 [--seed S]   draw the next batch, excluding every target already in targets.csv
            --genes A,B,C       or name the batch yourself (e.g. the most-used targets)
  worklist  [--batch N]         after fetch_papers.py --extract --backend claude: one row per block,
                                flagged subject / biomarker_only / no_mention; records the papers
  triage    [PMC ...]           section titles and flags only, per paper, to decide what to read
  show      ID|PMC ... [--max N] [--flag subject]   print block text
  skip      "reason" ID|PMC ... a block id, or a bare PMCID = every block of it not yet done
  todo                          blocks still to read, per target
  log                           rebuild evidence.csv, term_occurrences.csv, targets.csv counts
  summary                       write crawl/summaries/pass_NN.md for the periodic review pass

Nothing here edits data/. Findings reach the app only through a review chunk the owner applies.
"""
import argparse
import collections
import csv
import glob
import json
import os
import random
import re
import sys
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import literature_stats as ls  # noqa: E402
import readout_gaps as rg  # noqa: E402
import target_subject as ts  # noqa: E402

CRAWL = os.path.join(ROOT, "crawl")
PAPERS = os.path.join(ROOT, "build_output", "papers")
QUEUE = os.path.join(PAPERS, "agent_queue")
SURVEY = os.path.join(PAPERS, "survey.csv")

TARGETS = os.path.join(CRAWL, "targets.csv")
PAPER_LOG = os.path.join(CRAWL, "papers.csv")
SKIPPED = os.path.join(CRAWL, "skipped.csv")
EVIDENCE = os.path.join(CRAWL, "evidence.csv")
OCCURRENCES = os.path.join(CRAWL, "term_occurrences.csv")
CATALOG = os.path.join(CRAWL, "term_catalog.csv")
ALIASES = os.path.join(CRAWL, "alias_issues.csv")
WORKLISTS = os.path.join(CRAWL, "worklists")
SUMMARIES = os.path.join(CRAWL, "summaries")

TARGET_COLS = ["gene", "batch", "picked_on", "primary_class", "status", "papers", "papers_with_subject",
               "blocks", "records", "skipped", "new_entry_links", "new_readouts", "new_terms", "notes"]
PAPER_COLS = ["pmcid", "gene", "batch", "year", "title_verdict", "score", "n_blocks", "subject_blocks", "title"]
SKIP_COLS = ["assay_id", "reason", "read_on"]
EVIDENCE_COLS = ["gene", "batch", "pmcid", "assay_id", "confidence", "scope", "function_labels",
                 "technology_labels", "new_technology", "effect_direction", "pharmacology",
                 "protein_classes", "entries", "new_entry_links", "new_readouts", "quote", "title"]
OCC_COLS = ["catalog_term", "field", "verbatim", "nearest_controlled_value", "gene", "batch", "pmcid",
            "assay_id", "quote"]
# `variants` are paper verbatims and are matched EXACTLY by catalog_lookup(), so nothing phrased
# for the labeller belongs there. label_synonyms / label_suppresses carry the ChEMBL-side wording
# through to scripts/term_table.csv instead; see tools/papers_to_review.py.
CATALOG_COLS = ["term_id", "proposed_term", "kind", "nearest_controlled_value", "variants",
                "label_synonyms", "label_suppresses", "definition", "status", "first_pass",
                "notes"]
WORK_COLS = ["assay_id", "pmcid", "gene", "chars", "section", "flag", "block_path", "write_to"]
# One row per alias that cost papers. Written by `alias` as it is found, then enriched by `log`:
# the fix is applied in one pass over several batches, because changing an alias invalidates the
# cached Europe PMC search for every target already crawled.
ALIAS_COLS = ["gene", "alias", "kind", "collides_with", "papers_lost", "hits", "alias_only_hits",
              "alias_only_share", "found_in_batch", "action", "note"]
ALIAS_KINDS = ("alias", "species", "family", "mutation_notation")

# Pool: proteins whose papers exist and are about them. Prevalence does not matter here, so the
# campaign-count and readout-count filters of the yield sample are dropped.
VOCAB_FIELDS = {"technology_labels", "function_labels", "assay_category"}
MIN_HITS = 20
MIN_SUBJECT = 0.7
DEFAULT_SEED = 20261116


# ---------------------------------------------------------------- io

def read_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))


def write_rows(path, cols, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def append_rows(path, cols, rows):
    new = not os.path.exists(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=cols, extrasaction="ignore")
        if new:
            w.writeheader()
        for r in rows:
            w.writerow(r)


def num(x, cast=float):
    try:
        return cast(x)
    except (TypeError, ValueError):
        return cast(0)


def utf8_stdout():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


def survey_rows():
    return {r["gene"]: r for r in read_rows(SURVEY)}


def primary_class(row):
    return ((row.get("classes") or "").split(";")[0].strip()) or "unclassified"


def current_batch(targets):
    return max([num(t["batch"], int) for t in targets] or [0])


# ---------------------------------------------------------------- pick

def cmd_pick(args):
    targets = read_rows(TARGETS)
    done = {t["gene"] for t in targets}
    survey = survey_rows()
    batch = current_batch(targets) + 1
    today = datetime.now().strftime("%Y-%m-%d")
    if args.genes:
        # Named batch: the owner chose the targets (e.g. the most-used ones), so the pool
        # filters do not apply. Only "already searched" and "no survey row" are refused.
        picked, bad = [], []
        for g in [x.strip().upper() for x in args.genes.split(",") if x.strip()]:
            if g in done:
                bad.append("%s: already searched (batch %s)" %
                           (g, next(t["batch"] for t in targets if t["gene"] == g)))
            elif g not in survey:
                bad.append("%s: no survey row" % g)
            else:
                picked.append(survey[g])
        if bad:
            raise SystemExit("refusing to pick:\n   " + "\n   ".join(bad))
        append_rows(TARGETS, TARGET_COLS, [{"gene": r["gene"], "batch": batch, "picked_on": today,
                                            "primary_class": primary_class(r), "status": "picked"}
                                           for r in picked])
        genes = ",".join(sorted(r["gene"] for r in picked))
        print("batch %d: %d named targets" % (batch, len(picked)))
        for r in sorted(picked, key=lambda r: (primary_class(r), r["gene"])):
            print("   %-10s %-24s hits %s  subject %s" %
                  (r["gene"], primary_class(r)[:24], r.get("epmc_hits"), r.get("frac_subject")))
        print("\nnext:\n   py scripts\\fetch_papers.py --genes %s --extract --backend claude" % genes)
        print("   py scripts\\paper_crawl.py worklist")
        return
    pool = [r for r in survey.values()
            if r["gene"] not in done and not r.get("survey_error")
            and num(r.get("epmc_hits")) >= MIN_HITS and num(r.get("frac_subject")) >= MIN_SUBJECT]
    if not pool:
        raise SystemExit("pool is empty: every eligible target has been searched")
    rng = random.Random((args.seed or DEFAULT_SEED) + batch)
    by_class = collections.defaultdict(list)
    for r in pool:
        by_class[primary_class(r)].append(r)
    n = min(args.n, len(pool))
    # Stratified by primary class, largest remainder, so each batch mirrors the pool.
    quota = {c: n * len(v) / len(pool) for c, v in by_class.items()}
    alloc = {c: int(q) for c, q in quota.items()}
    for c in sorted(by_class, key=lambda c: -(quota[c] - alloc[c]))[:n - sum(alloc.values())]:
        alloc[c] += 1
    picked = []
    for c in sorted(by_class):
        if alloc[c]:
            picked += rng.sample(sorted(by_class[c], key=lambda r: r["gene"]), alloc[c])
    append_rows(TARGETS, TARGET_COLS, [{"gene": r["gene"], "batch": batch, "picked_on": today,
                                        "primary_class": primary_class(r), "status": "picked"}
                                       for r in picked])
    genes = ",".join(sorted(r["gene"] for r in picked))
    print("batch %d: %d targets from a pool of %d not yet searched" % (batch, len(picked), len(pool)))
    for r in sorted(picked, key=lambda r: (primary_class(r), r["gene"])):
        print("   %-10s %-24s hits %s" % (r["gene"], primary_class(r)[:24], r.get("epmc_hits")))
    print("\nnext:\n   py scripts\\fetch_papers.py --genes %s --extract --backend claude" % genes)
    print("   py scripts\\paper_crawl.py worklist")


# ---------------------------------------------------------------- worklist

def task_rows():
    path = os.path.join(QUEUE, "TASKS.md")
    out = []
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as h:
        for line in h:
            if not line.startswith("|") or "`" not in line:
                continue
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) < 7:
                continue
            _done, aid, gene, _chars, section, bpath, wto = cells[:7]
            out.append({"assay_id": aid.strip("`"), "gene": gene, "section": section,
                        "block_path": bpath.strip("`"), "write_to": wto.strip("`")})
    return out


def cmd_worklist(args):
    targets = read_rows(TARGETS)
    batch = args.batch if args.batch is not None else current_batch(targets)
    genes = {t["gene"] for t in targets if num(t["batch"], int) == batch}
    if not genes:
        raise SystemExit("no targets in batch %d" % batch)
    survey = survey_rows()
    pats = {}
    rows = []
    for t in task_rows():
        if t["gene"] not in genes:
            continue
        g = t["gene"]
        if g not in pats:
            al = ts.specific_aliases(rg.target_aliases(g, (survey.get(g) or {}).get("accession", "")),
                                     rg.normalise(g))
            pats[g] = ts.alias_pattern(al)
        try:
            text = open(os.path.join(ROOT, t["block_path"]), encoding="utf-8").read()
        except OSError:
            continue
        norm = ts.normalise(t["section"] + " " + text)
        pat = pats[g]
        flag = ("subject" if pat and ts.subject_positions(norm, pat)
                else "biomarker_only" if pat and pat.search(norm) else "no_mention")
        rows.append(dict(t, pmcid=t["assay_id"].split("-")[0], chars=len(text), flag=flag))
    by_paper = collections.defaultdict(list)
    for r in rows:
        by_paper[r["pmcid"]].append(r)
    # One paper can be queued under several genes of the same batch (a CYP panel, a kinase panel).
    # Each gene pass numbers that paper's blocks from A1, so the same id can name two different
    # blocks and both would be written to the same file. Give every distinct block one id, shared
    # by the genes that queued it; a block whose id is already taken gets the next free number.
    for p, rs in by_paper.items():
        taken, byid = {}, {}
        for r in sorted(rs, key=lambda r: (num(r["assay_id"].split("-A")[-1], int), r["gene"])):
            path = r["block_path"]
            if path in byid:
                r["assay_id"] = byid[path]
                continue
            aid = r["assay_id"]
            if aid in taken:
                n = max(num(x.split("-A")[-1], int) for x in taken) + 1
                aid = "%s-A%d" % (p, n)
            taken[aid] = path
            byid[path] = aid
            r["assay_id"] = aid
        for r in rs:
            r["write_to"] = os.path.join(QUEUE.replace(ROOT + os.sep, ""), p,
                                         "A%s.json" % r["assay_id"].split("-A")[-1])
    order = sorted(by_paper, key=lambda p: (by_paper[p][0]["gene"],
                                            -sum(r["flag"] == "subject" for r in by_paper[p])))
    out = []
    for p in order:
        out += sorted(by_paper[p], key=lambda r: num(r["assay_id"].split("-A")[-1], int))
    path = os.path.join(WORKLISTS, "batch_%02d.csv" % batch)
    write_rows(path, WORK_COLS, out)

    # Record the papers, because index.json is overwritten by every fetch.
    index = ls.read_json(os.path.join(PAPERS, "index.json"), {}) or {}
    known = {(r["pmcid"], r["gene"]) for r in read_rows(PAPER_LOG)}
    new = []
    for p in index.get("papers") or []:
        if p.get("gene") in genes and (p["pmcid"], p["gene"]) not in known:
            new.append({"pmcid": p["pmcid"], "gene": p["gene"], "batch": batch, "year": p.get("year"),
                        "title_verdict": p.get("title_verdict"), "score": p.get("score"),
                        "n_blocks": len(by_paper.get(p["pmcid"], [])),
                        "subject_blocks": sum(r["flag"] == "subject" for r in by_paper.get(p["pmcid"], [])),
                        "title": (p.get("title") or "")[:200]})
    append_rows(PAPER_LOG, PAPER_COLS, new)
    for t in targets:
        if t["gene"] in genes and t["status"] == "picked":
            t["status"] = "fetched"
    write_rows(TARGETS, TARGET_COLS, targets)
    c = collections.Counter(r["flag"] for r in out)
    print("batch %d worklist: %d blocks in %d papers  %s" % (batch, len(out), len(by_paper), dict(c)))
    print("   %s" % os.path.relpath(path, ROOT))
    print("targets with no queued paper: %s" % (", ".join(sorted(genes - {r["gene"] for r in out})) or "none"))


def load_worklist(batch=None):
    rows = []
    for path in sorted(glob.glob(os.path.join(WORKLISTS, "batch_*.csv"))):
        b = num(re.sub(r"\D", "", os.path.basename(path)), int)
        if batch is None or b == batch:
            rows += [dict(r, batch=b) for r in read_rows(path)]
    return rows


def skipped_ids():
    return {r["assay_id"] for r in read_rows(SKIPPED)}


def is_done(row, skipped):
    return row["assay_id"] in skipped or os.path.exists(os.path.join(ROOT, row["write_to"]))


# ---------------------------------------------------------------- reading helpers

def select(ids, rows):
    out = []
    for i in ids:
        out += [r for r in rows if r["assay_id"] == i or (re.fullmatch(r"PMC\d+", i) and r["pmcid"] == i)]
    return out


def cmd_triage(args):
    utf8_stdout()
    rows = load_worklist()
    sk = skipped_ids()
    papers = {r["pmcid"]: r for r in read_rows(PAPER_LOG)}
    pend = [r for r in rows if not is_done(r, sk)]
    if args.ids:
        pend = select(args.ids, pend)
    last = None
    for r in pend:
        if r["pmcid"] != last:
            last = r["pmcid"]
            print("\n## %s  %s  %s" % (r["pmcid"], r["gene"], (papers.get(r["pmcid"], {}).get("title") or "")[:110]))
        print("   %-16s %-14s %5s  %s" % (r["assay_id"], r["flag"], r["chars"], r["section"][:80]))


def cmd_show(args):
    utf8_stdout()
    rows = load_worklist()
    sk = skipped_ids()
    for r in select(args.ids, rows):
        if is_done(r, sk) and not args.all:
            continue
        if args.flag and r["flag"] != args.flag:
            continue
        txt = " ".join(open(os.path.join(ROOT, r["block_path"]), encoding="utf-8").read().split())
        print("\n--- %s %s [%s] %s (%d chars) -> %s" % (r["assay_id"], r["gene"], r["flag"], r["section"],
                                                       len(txt), r["write_to"]))
        print(txt[:args.max] + (" [...truncated]" if len(txt) > args.max else ""))


def cmd_skip(args):
    rows = load_worklist()
    sk = skipped_ids()
    ids = [r["assay_id"] for r in select(args.ids, rows) if not is_done(r, sk)]
    ids += [i for i in args.ids if not re.fullmatch(r"PMC\d+", i) and i not in {r["assay_id"] for r in rows}]
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    append_rows(SKIPPED, SKIP_COLS, [{"assay_id": i, "reason": args.reason, "read_on": now} for i in ids])
    print("%d block(s) skipped: %s" % (len(ids), args.reason))


def alias_share(gene, alias, survey):
    """How much of this target's search rests on this one alias: hits that name the alias but
    neither the gene symbol nor the full protein name, over all hits.

    Measured from the cached search payload, so it costs nothing and needs no full text. A share
    near 1 means the alias is carrying the whole search and the papers are somebody else's subject
    (UNC119 'rg4' = 0.92). A middling share does NOT settle it: a legitimate common name scores the
    same as a real collision (GPBAR1 'tgr5' = 0.48, F13A1 'f13a' = 0.42), which is why the verdict
    is recorded per target rather than derived from a threshold."""
    row = survey.get(gene) or {}
    path = os.path.join(ROOT, "build_cache", "papers", "search", "%s.json" % gene)
    if not os.path.exists(path):
        return None
    al = ts.specific_aliases(rg.target_aliases(gene, row.get("accession", "")), rg.normalise(gene))
    a = ts.normalise(alias)
    if a == ts.normalise(gene):
        return None  # the share compares an alias against the symbol; for the symbol it is vacuous
    rest = [x for x in al if x != a]
    pr, pa = ts.alias_pattern(rest), ts.alias_pattern([a])
    if not pa:
        return None
    only = n = 0
    for r in (ls.read_json(path, {}) or {}).get("results") or []:
        text = ts.normalise("%s %s" % (r.get("title") or "", r.get("abstractText") or ""))
        n += 1
        if pa.search(text) and not (pr and pr.search(text)):
            only += 1
    return (n, only, round(only / n, 2)) if n else None


def cmd_alias(args):
    rows = read_rows(ALIASES)
    targets = {t["gene"]: t for t in read_rows(TARGETS)}
    key = (args.gene, ts.normalise(args.alias))
    row = next((r for r in rows if (r["gene"], ts.normalise(r["alias"])) == key), None)
    if row is None:
        row = {"gene": args.gene, "alias": args.alias,
               "found_in_batch": (targets.get(args.gene) or {}).get("batch", ""), "action": "open"}
        rows.append(row)
    row["kind"] = args.kind
    row["collides_with"] = args.collides_with or row.get("collides_with", "")
    if args.note:
        row["note"] = args.note
    write_rows(ALIASES, ALIAS_COLS, rows)
    print("logged %s alias '%s' (%s)" % (args.gene, args.alias, args.kind))


def enrich_aliases(survey, skips):
    """Fill the measured columns and the papers-lost count. Called by `log` so the file stays
    current as more batches land, without re-deciding any verdict."""
    rows = read_rows(ALIASES)
    if not rows:
        return rows
    wl = {r["assay_id"]: r for r in load_worklist()}
    lost = collections.defaultdict(set)
    for s in skips:
        if "collision" in (s.get("reason") or "").lower() and s["assay_id"] in wl:
            lost[wl[s["assay_id"]]["gene"]].add(wl[s["assay_id"]]["pmcid"])
    for r in rows:
        r["papers_lost"] = len(lost.get(r["gene"], ()))
        got = alias_share(r["gene"], r["alias"], survey)
        if got:
            r["hits"], r["alias_only_hits"], r["alias_only_share"] = got
    write_rows(ALIASES, ALIAS_COLS, rows)
    return rows


def squash(text):
    return re.sub(r"\s+", " ", text or "").strip()


def cmd_record(args):
    """Write lean assay records from one JSON list to their worklist paths. Each quote must occur
    verbatim (whitespace aside) in its block, or nothing is written for that record."""
    recs = json.load(open(args.spec, encoding="utf-8"))
    rows = {r["assay_id"]: r for r in load_worklist()}
    ok = bad = 0
    for a in recs:
        # A block that describes several assays: give the extra ones ids A101, A102... and name the
        # block they come from in "block"; they are written beside it and merged like any other.
        parent = a.pop("block", None)
        r = rows.get(parent or a.get("assay_id"))
        if not r:
            print("   %s: not in any worklist" % a.get("assay_id")); bad += 1; continue
        if parent:
            n = a["assay_id"].split("-A")[-1]
            r = dict(r, write_to=os.path.join(os.path.dirname(r["write_to"]), "A%s.json" % n))
        block = squash(open(os.path.join(ROOT, r["block_path"]), encoding="utf-8").read())
        a.setdefault("source_section", r["section"])
        a.setdefault("technology_labels", [])
        a.setdefault("unmapped_terms", [])
        a.setdefault("extraction_flags", [])
        for e in a.get("evidence") or []:
            e.setdefault("section", r["section"])
        quotes = [e.get("quote") for e in a.get("evidence") or []] + [u.get("quote") for u in a["unmapped_terms"]]
        missing = [q for q in quotes if q and squash(q) not in block]
        if missing or not a.get("evidence"):
            print("   %s: quote not found in block: %s" % (a["assay_id"], (missing or ["(no evidence)"])[0][:90]))
            bad += 1
            continue
        out = os.path.join(ROOT, r["write_to"])
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w", encoding="utf-8") as h:
            json.dump(a, h, indent=2, ensure_ascii=False)
        ok += 1
    print("%d record(s) written, %d rejected" % (ok, bad))


def cmd_todo(args):
    rows = load_worklist()
    sk = skipped_ids()
    left = collections.Counter()
    total = collections.Counter()
    for r in rows:
        total[(r["batch"], r["gene"])] += 1
        if not is_done(r, sk):
            left[(r["batch"], r["gene"])] += 1
    for k in sorted(total):
        if left[k]:
            print("batch %d  %-10s %3d of %3d blocks left" % (k[0], k[1], left[k], total[k]))
    print("total left: %d" % sum(left.values()))


# ---------------------------------------------------------------- log

def generic_tech_ok(entry, techs):
    """Generic entries (G1-G7) are defined by one technique: a binding result by MST must not link
    the DSF, SPR, ITC and MS entries just because they share 'Ligand binding affinity'."""
    if "generic" not in entry["classes"]:
        return True
    # load_entries lower-cases the readouts; lower again so a caller need not know that.
    pieces = [p.strip() for p in re.split(r"[;,]", " ".join(entry["readouts"]).lower()) if p.strip()]
    for t in techs:
        tl = t.lower()
        head = re.split(r"[ /(]", tl)[0]
        for p in pieces:
            if p in tl or tl in p or (len(head) > 2 and head in p):
                return True
    return False


def catalog_lookup():
    look = {}
    for c in read_rows(CATALOG):
        for v in (c.get("variants") or "").split("|"):
            if v.strip():
                look[v.strip().lower()] = c["proposed_term"]
    return look


def linked_today():
    out = collections.defaultdict(set)
    for u in read_rows(os.path.join(ROOT, "data", "usage_stats.csv")):
        if u["scope_type"] == "gene" and not u["readout"]:
            out[u["scope_id"]].add(u["entry_id"])
    return out


def cmd_log(args):
    survey = survey_rows()
    entries = ls.load_entries()
    camp = rg.campaign_data()
    linked = linked_today()
    batch_of = {}
    for p in read_rows(PAPER_LOG):
        batch_of.setdefault(p["pmcid"], p["batch"])
    titles = {p["pmcid"]: p["title"] for p in read_rows(PAPER_LOG)}
    look = catalog_lookup()
    ev, occ = [], []
    for path in sorted(glob.glob(os.path.join(PAPERS, "PMC*.json"))):
        d = ls.read_json(path, {}) or {}
        pm = (d.get("paper") or {}).get("pmcid")
        if pm not in batch_of:
            continue
        b = batch_of[pm]
        paper_gene = ((d.get("paper") or {}).get("target_genes") or [""])[0]
        # --agent-merge hoists each assay's unmapped_terms to the record, keyed by assay_id
        unm_by = collections.defaultdict(list)
        for u in d.get("unmapped_terms") or []:
            if u.get("field") in VOCAB_FIELDS:
                unm_by[u.get("assay_id")].append(u.get("verbatim") or "")
        for a in d.get("assays") or []:
            g = (a.get("target") or {}).get("gene_symbol") or paper_gene
            techs = a.get("technology_labels") or []
            cell = a.get("scope") == "cell-based"
            classes = ls.classes_of(g, survey)
            ents = [e["entry_id"] for e in ls.entries_for(entries, classes, a.get("function_labels") or [], cell)
                    if generic_tech_ok(e, techs)]
            known = set((camp.get(g) or {}).get("readouts") or ())
            unm = unm_by.get(a.get("assay_id"), []) + [u.get("verbatim", "") for u in a.get("unmapped_terms") or []
                                                         if u.get("field") in VOCAB_FIELDS]
            ev.append({"gene": g, "batch": b, "pmcid": pm, "assay_id": a.get("assay_id"),
                       "confidence": a.get("confidence"), "scope": a.get("scope"),
                       "function_labels": "; ".join(a.get("function_labels") or []),
                       "technology_labels": "; ".join(techs),
                       "new_technology": "; ".join(look.get(v.lower(), v) for v in unm),
                       "effect_direction": a.get("effect_direction") or "",
                       "pharmacology": "; ".join(p for p in a.get("pharmacology") or [] if p != "unclear"),
                       "protein_classes": "; ".join(classes) or "none",
                       "entries": "; ".join(ents),
                       "new_entry_links": "; ".join(e for e in ents if e not in linked[g]),
                       "new_readouts": "; ".join(t for t in techs if t not in known
                                                 and t not in rg.GENERIC_TECHNIQUES),
                       "quote": ls.best_quote(a)[:400], "title": (titles.get(pm) or "")[:120]})
        for u in d.get("unmapped_terms") or []:
            if u.get("field") not in VOCAB_FIELDS:
                continue  # schema-detail notes (buffer, host, derived DMSO), not vocabulary
            aid = u.get("assay_id")
            a = next((x for x in d.get("assays") or [] if x.get("assay_id") == aid), {})
            g = (a.get("target") or {}).get("gene_symbol") or paper_gene
            v = (u.get("verbatim") or "").strip()
            occ.append({"catalog_term": look.get(v.lower(), ""), "field": u.get("field"), "verbatim": v,
                        "nearest_controlled_value": u.get("nearest_controlled_value") or "",
                        "gene": g, "batch": b, "pmcid": pm, "assay_id": aid,
                        "quote": " ".join((u.get("quote") or "").split())[:300]})
    write_rows(EVIDENCE, EVIDENCE_COLS, ev)
    write_rows(OCCURRENCES, OCC_COLS, occ)
    if not os.path.exists(CATALOG):
        write_rows(CATALOG, CATALOG_COLS, [])

    # target counts and status
    targets = read_rows(TARGETS)
    wl = load_worklist()
    sk = skipped_ids()
    papers = read_rows(PAPER_LOG)
    hm = [e for e in ev if e["confidence"] in ("high", "medium")]
    for t in targets:
        g = t["gene"]
        mine = [r for r in wl if r["gene"] == g]
        if not mine and t["status"] in ("picked",):
            continue
        t["papers"] = sum(1 for p in papers if p["gene"] == g)
        t["papers_with_subject"] = sum(1 for p in papers if p["gene"] == g and num(p["subject_blocks"], int) > 0)
        t["blocks"] = len(mine)
        t["skipped"] = sum(1 for r in mine if r["assay_id"] in sk)
        t["records"] = sum(1 for e in ev if e["gene"] == g)
        t["new_entry_links"] = "; ".join(sorted({x for e in hm if e["gene"] == g
                                                 for x in e["new_entry_links"].split("; ") if x}))
        t["new_readouts"] = "; ".join(sorted({x for e in hm if e["gene"] == g
                                              for x in e["new_readouts"].split("; ") if x}))
        t["new_terms"] = len({o["verbatim"].lower() for o in occ if o["gene"] == g})
        if mine and all(is_done(r, sk) for r in mine):
            t["status"] = "read"
        elif mine:
            t["status"] = "reading"
        elif t["status"] == "fetched":
            t["status"] = "no_papers"
    write_rows(TARGETS, TARGET_COLS, targets)
    al = enrich_aliases(survey, read_rows(SKIPPED))
    if al:
        print("alias issues: %d logged, %d still open, %d paper(s) lost to collisions"
              % (len(al), sum(1 for r in al if r.get("action") == "open"),
                 sum(num(r.get("papers_lost"), int) for r in al)))
    st = collections.Counter(t["status"] for t in targets)
    print("evidence rows %d (%d high/medium), term occurrences %d (%d not yet catalogued)"
          % (len(ev), len(hm), len(occ), sum(1 for o in occ if not o["catalog_term"])))
    print("targets: %s" % ", ".join("%s %d" % kv for kv in sorted(st.items())))


# ---------------------------------------------------------------- summary (periodic pass)

def cmd_summary(args):
    os.makedirs(SUMMARIES, exist_ok=True)
    prev = sorted(glob.glob(os.path.join(SUMMARIES, "pass_*.md")))
    n = len(prev) + 1
    last_batch = -1  # batch 0 is the yield sample; the first pass covers it
    if prev:
        m = re.search(r"batches_through:\s*(\d+)", open(prev[-1], encoding="utf-8").read())
        last_batch = num(m.group(1), int) if m else 0
    targets = read_rows(TARGETS)
    through = max([num(t["batch"], int) for t in targets if t["status"] == "read"] or [0])
    if through <= last_batch:
        raise SystemExit("batch %d is already covered by %s; read the next batch first"
                         % (through, os.path.basename(prev[-1])))
    ev = read_rows(EVIDENCE)
    occ = read_rows(OCCURRENCES)
    cat = read_rows(CATALOG)

    L = ["# Paper crawl - pass %d" % n, "",
         "<!-- batches_through: %d -->" % through,
         "Written %s. Covers batches %d-%d. Targets read so far: %d of %d picked." % (
             datetime.now().strftime("%Y-%m-%d"), last_batch + 1, through,
             sum(t["status"] == "read" for t in targets), len(targets)), ""]

    # 1. vocabulary
    by_term = collections.defaultdict(lambda: {"papers": set(), "genes": set()})
    for o in occ:
        k = o["catalog_term"] or ("UNCATALOGUED: " + o["verbatim"])
        by_term[k]["papers"].add(o["pmcid"])
        by_term[k]["genes"].add(o["gene"])
    L += ["## 1. Technologies and assay types not in the vocabulary", "",
          "| proposed term | kind | status | papers | targets | targets named |", "|---|---|---|---|---|---|"]
    for c in sorted(cat, key=lambda c: -len(by_term[c["proposed_term"]]["papers"])):
        s = by_term[c["proposed_term"]]
        L.append("| %s | %s | %s | %d | %d | %s |" % (c["proposed_term"], c["kind"], c["status"], len(s["papers"]),
                                                     len(s["genes"]), ", ".join(sorted(s["genes"]))[:80]))
    unc = sorted(k for k in by_term if k.startswith("UNCATALOGUED: "))
    if unc:
        L += ["", "Not yet catalogued (assign each to a catalog term, or add a new term):", ""]
        L += ["- %s (%s)" % (k[14:], ", ".join(sorted(by_term[k]["genes"]))) for k in unc]

    # 2. target-assay evidence from this pass
    new = [e for e in ev if num(e["batch"], int) > last_batch and e["confidence"] in ("high", "medium")]
    L += ["", "## 2. New evidence that an assay type is used for a target (batches %d-%d)" % (last_batch + 1, through), "",
          "Only rows the app does not already show: an entry the target is not linked to today, a readout "
          "no campaign states, or a technology with no term yet.", "",
          "| target | class | new entry link | new readout | new technology | direction | paper |",
          "|---|---|---|---|---|---|---|"]
    k = 0
    for e in sorted(new, key=lambda e: (e["gene"], e["pmcid"])):
        if not (e["new_entry_links"] or e["new_readouts"] or e["new_technology"]):
            continue
        k += 1
        L.append("| %s | %s | %s | %s | %s | %s | %s |" % (e["gene"], e["protein_classes"], e["new_entry_links"],
                                                          e["new_readouts"], e["new_technology"][:60],
                                                          e["effect_direction"], e["pmcid"]))
    nocls = sorted({e["gene"] for e in new if e["protein_classes"] == "none"})
    L += ["", "%d rows. Targets with evidence but no protein class (no entry can link to them): %s"
          % (k, ", ".join(nocls) or "none"), ""]

    # 3. per-target outcome of this pass
    L += ["## 3. Targets read in this pass", "", "| target | class | papers | records | new links | new readouts | new terms |",
          "|---|---|---|---|---|---|---|"]
    for t in sorted(targets, key=lambda t: t["gene"]):
        if t["status"] == "read" and last_batch < num(t["batch"], int) <= through:
            L.append("| %s | %s | %s | %s | %s | %s | %s |" % (t["gene"], t["primary_class"], t["papers"], t["records"],
                                                              t["new_entry_links"], t["new_readouts"], t["new_terms"]))
    L += ["", "## 4. Observations", "", "(written by the agent: patterns, systemic issues, terms ready to propose)", ""]
    path = os.path.join(SUMMARIES, "pass_%02d.md" % n)
    with open(path, "w", encoding="utf-8") as h:
        h.write("\n".join(L))
    print("wrote %s  (%d uncatalogued term(s) to assign)" % (os.path.relpath(path, ROOT), len(unc)))


# ---------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pick"); p.add_argument("--n", type=int, default=20); p.add_argument("--seed", type=int)
    p.add_argument("--genes", help="comma-separated gene symbols: pick these instead of a random draw")
    p = sub.add_parser("worklist"); p.add_argument("--batch", type=int)
    p = sub.add_parser("triage"); p.add_argument("ids", nargs="*")
    p = sub.add_parser("show"); p.add_argument("ids", nargs="+"); p.add_argument("--max", type=int, default=3000)
    p.add_argument("--flag"); p.add_argument("--all", action="store_true", help="include blocks already done")
    p = sub.add_parser("skip"); p.add_argument("reason"); p.add_argument("ids", nargs="+")
    p = sub.add_parser("record"); p.add_argument("spec", help="JSON file: a list of lean assay objects")
    p = sub.add_parser("alias"); p.add_argument("gene"); p.add_argument("alias")
    p.add_argument("--kind", choices=ALIAS_KINDS, required=True)
    p.add_argument("--collides-with", help="what the alias actually matched")
    p.add_argument("--note")
    sub.add_parser("todo")
    sub.add_parser("log")
    sub.add_parser("summary")
    args = ap.parse_args(argv)
    {"pick": cmd_pick, "worklist": cmd_worklist, "triage": cmd_triage, "show": cmd_show, "skip": cmd_skip,
     "todo": cmd_todo, "log": cmd_log, "record": cmd_record, "alias": cmd_alias, "summary": cmd_summary}[args.cmd](args)


if __name__ == "__main__":
    main()
