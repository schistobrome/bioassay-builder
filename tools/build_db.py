#!/usr/bin/env python3
"""Check the data/ tables and build docs/assay_db.json for the website.

Run from the repo folder:
    python tools/build_db.py                 # check + build docs/assay_db.json
    python tools/build_db.py --include-examples   # also load examples/usage_stats.example.csv (demo only)
Standard library only. Problems are listed as ERROR (build stops) or WARNING (build continues).
"""
import argparse, csv, json, os, re, sys
from datetime import date

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
CATEGORIES = {"biochemical_hts", "cell_functional", "biophysical_binding", "orthogonal_confirmation", "mechanistic_counterscreen"}
ROLES = {"primary", "confirmation", "counterscreen"}
LEVELS = {"very_low", "low", "moderate", "high", "very_high"}
STAGES = {"primary_screen", "confirmation", "counterscreen", "any"}
# inhibition_labeled / activation_with_unspecified: the same directions with unspecified campaigns left out / added
USAGE_EFFECTS = {"any", "inhibition", "activation", "inhibition_labeled", "activation_with_unspecified"}
USAGE_COLUMNS = ["entry_id", "scope_type", "scope_id", "stage", "effect", "readout", "plate", "total", "classifiable",
                 "using", "lit_total", "lit_using", "score", "level", "low_confidence", "source"]
ASSAY_EFFECTS = {"", "both", "inhibition", "activation"}
SCOPES = {"cell-free", "cell-based"}
# Curated from the NIH Assay Guidance Manual (review chunk 07); the website holds the wording for each code.
QC_METRICS = {"", "hts", "hcs", "potency"}
INTERFERENCE = {"luciferase", "fluorescence", "absorbance", "proximity", "aggregation", "thiol_redox", "cytotoxicity",
                "imaging", "label_free"}
MAX_JSON_MB = 15          # docs/assay_db.json above this size gets a warning
# The labeller's vocabulary. Its technology terms are the readout labels usage_stats.csv can carry.
TERM_TABLE = os.path.join(ROOT, "scripts", "term_table.csv")
# Regex syntax Python accepts but browsers do not. readout_labels.csv patterns run in the page as
# new RegExp(pattern, "i"), so one of these would compile here and silently fail there.
PY_ONLY_REGEX = [(r"\(\?P", "(?P...) named groups (write (?<name>...))"),
                 (r"\(\?#", "(?#...) comments"),
                 (r"\(\?[aiLmsux]+[):]", "inline flags such as (?i) (the page already ignores case)"),
                 (r"\(\?>", "atomic groups (?>...)"),
                 (r"[*+?}][+]", "possessive quantifiers such as *+"),
                 (r"\\[AZ]", r"\A or \Z (use ^ and $)")]
ERRORS, WARNINGS = [], []


def err(msg): ERRORS.append(msg)
def warn(msg): WARNINGS.append(msg)


def check_readout_labels(rows, usage_counts, term_tech):
    """Check data/readout_labels.csv against the data and the vocabulary.

    The page shows a readout's campaign share next to an assay's curated readouts by finding which
    labels the curated text names: each row here is a label and the words (a regex) that name it.
    Three things must agree, and before this table existed nothing made them:
      usage_counts  {label: n rows} -- readout labels usage_stats.csv actually carries
      term_tech     set of technology terms in scripts/term_table.csv, or None if it is absent
      rows          this table
    A label in the data with no row here cannot be matched at all: ERROR. A term in the vocabulary
    with no row here is harmless until the relabel puts it in the data, so it is a WARNING that
    says what to add before then.
    Returns (errors, warnings, [[label, pattern], ...] in file order)."""
    errors, warnings, out, seen = [], [], [], set()
    for i, r in enumerate(rows, start=2):
        label, pat = (r.get("label") or "").strip(), (r.get("pattern") or "").strip()
        where = "readout_labels.csv line %d" % i
        if not label:
            errors.append("%s: no label" % where); continue
        where = "readout_labels.csv '%s'" % label
        if label in seen:
            errors.append("%s: duplicate label" % where); continue
        seen.add(label)
        if not pat:
            errors.append("%s: no pattern, so no curated readout can ever match it" % where); continue
        try:
            re.compile(pat, re.I)
        except re.error as e:
            errors.append("%s: pattern does not compile (%s)" % (where, e)); continue
        bad = [why for rx, why in PY_ONLY_REGEX if re.search(rx, pat)]
        if bad:
            errors.append("%s: pattern uses %s, which works here but not in the browser"
                          % (where, "; ".join(bad))); continue
        out.append([label, pat])
    for label in sorted(set(usage_counts) - seen):
        errors.append("usage_stats.csv carries readout label '%s' (%d rows) but data/readout_labels.csv "
                      "has no row for it, so the website cannot match it to any assay's curated readouts. "
                      "Add: %s,<regex for the words assays.csv uses for it>,"
                      % (label, usage_counts[label], label))
    if term_tech is not None:
        for t in sorted(term_tech - seen):
            if t not in usage_counts:      # already an ERROR above when it is in the data
                warnings.append("scripts/term_table.csv has technology term '%s' but data/readout_labels.csv "
                                "has no row for it. Add one before the m2 relabel, or this build will fail "
                                "once the usage data carries that label." % t)
        for t in sorted(seen - term_tech):
            warnings.append("readout_labels.csv '%s' is not a technology term in scripts/term_table.csv, "
                            "so no campaign can carry it (renamed or removed there?)" % t)
    return errors, warnings, out


def read(name, path=None):
    path = path or os.path.join(DATA, name)
    with open(path, newline="", encoding="utf-8-sig") as f:
        return [{k: (v or "").strip() for k, v in r.items() if k} for r in csv.DictReader(f)]


def lst(cell):
    return [x.strip() for x in re.split(r"[;|]", cell or "") if x.strip()]


def num(cell, where):
    if cell == "":
        return None
    try:
        return int(float(cell))
    except ValueError:
        err("%s: '%s' is not a number" % (where, cell))
        return None


def level_for(score, levels):
    for key, floor in levels:
        if score >= floor:
            return key
    return levels[-1][0]


def priority(text):
    nums = [int(x) for x in re.findall(r"P\s*(\d)", text or "")]
    return {"best": min(nums) if nums else None, "text": text}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-examples", action="store_true", help="demo build: add illustrative usage numbers")
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "assay_db.json"))
    args = ap.parse_args()

    config = json.load(open(os.path.join(DATA, "config.json"), encoding="utf-8"))
    classes = read("function_classes.csv")
    assays = read("assays.csv")
    kits = read("kits.csv")
    cites = read("citations.csv")
    usage = read("usage_stats.csv")
    if args.include_examples:
        usage += read("", os.path.join(ROOT, "examples", "usage_stats.example.csv"))

    usage_labels = {}
    for u in usage:
        rd = (u.get("readout") or "").strip()
        if rd and not rd.startswith("plate:"):      # plate-format rows share the readout column
            usage_labels[rd] = usage_labels.get(rd, 0) + 1
    term_tech = None
    if os.path.exists(TERM_TABLE):
        term_tech = {r["term"] for r in read("", TERM_TABLE) if r.get("axis", "").lower() == "technology"}
    else:
        warn("scripts/term_table.csv not found, so readout labels are only checked against the usage data")
    e_rl, w_rl, readout_labels = check_readout_labels(read("readout_labels.csv"), usage_labels, term_tech)
    ERRORS.extend(e_rl); WARNINGS.extend(w_rl)

    cond_path = os.path.join(DATA, "condition_stats.csv")        # optional, written by build_database.py
    conds = []
    known_entries = {a["entry_id"] for a in assays}             # built once (was rebuilt for every condition row)
    if os.path.exists(cond_path):
        for u in read("condition_stats.csv"):
            if u["entry_id"] not in known_entries:
                continue
            def fnum(x): return float(x) if x not in ("", None) else None
            conds.append({"entry_id": u["entry_id"], "scope": {"type": u["scope_type"], "id": u["scope_id"]},
                          "parameter": u["parameter"], "kind": u["kind"], "n_using": num(u["n_campaigns_using"], "condition_stats.csv"),
                          "n_stating": num(u["n_campaigns_stating"], "condition_stats.csv"),
                          "median": fnum(u["median"]), "p25": fnum(u["p25"]), "p75": fnum(u["p75"]), "min": fnum(u["min"]), "max": fnum(u["max"]),
                          "top": [[v.rsplit("|", 1)[0], float(v.rsplit("|", 1)[1])] for v in [x.strip() for x in u["top_values"].split(";")] if "|" in v]})

    class_ids = {c["class_id"] for c in classes}
    if len(class_ids) != len(classes):
        err("function_classes.csv: duplicate class_id")
    entry_ids = set()
    for a in assays:
        e = a["entry_id"]
        if not e: err("assays.csv: a row has no entry_id"); continue
        if e in entry_ids: err("assays.csv: duplicate entry_id %s" % e)
        entry_ids.add(e)
        for c in lst(a["classes"]):
            if c not in class_ids: err("assays.csv %s: unknown class '%s'" % (e, c))
        if a["assay_category"] not in CATEGORIES: err("assays.csv %s: bad assay_category '%s'" % (e, a["assay_category"]))
        if a["scope"] not in SCOPES: err("assays.csv %s: scope must be cell-free or cell-based, got '%s'" % (e, a["scope"]))
        if not (lst(a["role"].replace(" or ", ";")) and lst(a["role"].replace(" or ", ";"))[0].lower() in ROLES):
            err("assays.csv %s: bad role '%s'" % (e, a["role"]))
        for col in ("screening_suitability", "biological_relevance"):
            if a[col] not in LEVELS: err("assays.csv %s: %s must be one of %s, got '%s'" % (e, col, sorted(LEVELS), a[col]))
        for d in lst(a.get("sources_doi", "")):
            if not re.match(r"^(10\.\d{4,9}/\S+|https?://\S+)$", d): warn("assays.csv %s: '%s' does not look like a DOI" % (e, d))
        if not a["status"]: warn("assays.csv %s: no review status" % e)
        if a.get("effect", "").lower() not in ASSAY_EFFECTS: err("assays.csv %s: effect must be blank, both, inhibition or activation, got '%s'" % (e, a["effect"]))
        if a.get("qc_metric", "") not in QC_METRICS: err("assays.csv %s: qc_metric must be blank or one of %s, got '%s'" % (e, sorted(QC_METRICS - {""}), a["qc_metric"]))
        for t in lst(a.get("interference", "")):
            if t not in INTERFERENCE: err("assays.csv %s: unknown interference code '%s' (allowed: %s)" % (e, t, ", ".join(sorted(INTERFERENCE))))
        if a.get("dmso_pct") and not re.match(r"^\d+(\.\d+)?(-\d+(\.\d+)?)?$", a["dmso_pct"]):
            err("assays.csv %s: dmso_pct must be a number or a range such as 1-2, got '%s'" % (e, a["dmso_pct"]))

    for k in kits:
        for e in lst(k["applies_to_entries"]):
            if e not in entry_ids: err("kits.csv %s: unknown entry '%s'" % (k["kit_id"], e))
        if not k["url"] and k["status"] != "unverified": warn("kits.csv %s: no url" % k["kit_id"])

    cite_by_doi = {c["doi"].lower(): c for c in cites if c["doi"]}
    w = config["prevalence"]["weights"]
    levels = config["prevalence"]["levels"]
    min_n = config["prevalence"]["min_campaigns_for_confidence"]
    # Usage rows are written as a compact table (schema 0.2): one short list per row, and the long, repeated
    # source / query / notes text stored once in a lookup list. docs/index.html turns them back into objects.
    stats, sources, source_index = [], [], {}
    for u in usage:
        where = "usage_stats.csv %s" % (u["stat_id"] or "?")
        if u["entry_id"] not in entry_ids: err("%s: unknown entry '%s'" % (where, u["entry_id"])); continue
        if u["scope_type"] not in ("gene", "function_class"): err("%s: scope_type must be gene or function_class" % where)
        if u["scope_type"] == "function_class" and u["scope_id"] not in class_ids: err("%s: unknown class '%s'" % (where, u["scope_id"]))
        if u["stage"] not in STAGES: err("%s: stage must be one of %s" % (where, sorted(STAGES)))
        if (u.get("effect") or "any") not in USAGE_EFFECTS: err("%s: effect must be one of %s" % (where, sorted(USAGE_EFFECTS)))
        nt, nc, nu = (num(u[c], where) for c in ("n_campaigns_total", "n_campaigns_classifiable", "n_campaigns_using"))
        lt, lu = num(u["n_literature_total"], where), num(u["n_literature_using"], where)
        if nc is not None and nu is not None and nu > nc: err("%s: n_campaigns_using is larger than n_campaigns_classifiable" % where)
        if lt is not None and lu is not None and lu > lt: err("%s: n_literature_using is larger than n_literature_total" % where)
        sc = nu / nc if nc else None
        sl = lu / lt if lt else None
        parts = [(w["campaigns"], sc), (w["literature"], sl)]
        parts = [(wt, s) for wt, s in parts if s is not None]
        score = sum(wt * s for wt, s in parts) / sum(wt for wt, s in parts) if parts else None
        rd = u["readout"]
        plate = rd[6:] if rd.startswith("plate:") else None      # plate-format rows share the readout column
        if plate is not None:
            rd = ""
        src = (u["source"], u["query"], u["retrieved_on"] or None, u["method_version"] or None, u["notes"])
        if src not in source_index:
            source_index[src] = len(sources)
            sources.append(list(src))
        stats.append([
            u["entry_id"], u["scope_type"], u["scope_id"], u["stage"], u.get("effect") or "any", rd or None, plate,
            nt, nc, nu, lt, lu,
            round(score, 4) if score is not None else None,
            level_for(score, levels) if score is not None else None,
            1 if (nc is None or nc < min_n) else 0,
            source_index[src]])

    out = {
        # 0.3: readout_labels added (was hard-coded as TECH in docs/index.html)
        "meta": {"schema_version": "0.3", "built_on": date.today().isoformat(),
                 "last_curated": max([a.get("reviewed_on", "") for a in assays if a.get("reviewed_on")] or [None]),
                 "counts": {"assays": len(assays), "kits": len(kits), "usage_stats": len(stats), "condition_stats": len(conds)},
                 "includes_example_data": bool(args.include_examples)},
        "config": config,
        "function_classes": [{
            "class_id": c["class_id"], "label": c["label"], "ec_prefixes": lst(c["ec_prefixes"]), "keywords": lst(c["keywords"]),
            "go_terms": lst(c["go_terms"]), "domain_patterns": lst(c["domain_patterns"]), "needs_domain": c["needs_domain"].lower() == "yes",
            "modality": c["modality"] or None, "excludes": lst(c["excludes"]), "notes": c["notes"]} for c in classes],
        "assays": [{
            "entry_id": a["entry_id"], "assay_name": a["assay_name"], "classes": lst(a["classes"]),
            "assay_category": a["assay_category"], "scope": a["scope"], "role": a["role"],
            "effect": (a.get("effect") or "").lower() or None,
            "development_priority": priority(a["development_priority"]),
            "screening_suitability": a["screening_suitability"], "biological_relevance": a["biological_relevance"],
            "rating_notes": a["rating_notes"], "measures": a["measures"], "readouts": lst(a["readouts"]),
            "evidence_labels": lst(a.get("evidence_labels", "")), "conditions": a["conditions"], "reagents": a["reagents"],
            # The biological preparation an assay runs in (a Ba/F3 line, HepaRG, a donor microsome
            # panel, an organoid). Separate from `scope`, which only says cell-free or cell-based,
            # and from `reagents`, which lists what goes into the well. Blank for most entries:
            # campaign descriptions rarely state it, so this fills from the literature.
            "systems": lst(a.get("system", "")),
            "rationale": a["rationale"], "plate_formats": lst(a.get("plate_formats", "")),
            "mechanism": {"learn_more_url": a.get("learn_more_url") or None},
            "quality": {"qc_metric": a.get("qc_metric") or None, "dmso_pct": a.get("dmso_pct") or None,
                        "interference": lst(a.get("interference", ""))},
            "sources": [{"doi": d, **({k: cite_by_doi[d.lower()][k] for k in ("title", "year", "journal") if cite_by_doi.get(d.lower())})} for d in lst(a.get("sources_doi", ""))],
            "review": {"status": a["status"], "reviewer": a.get("reviewer") or None, "reviewed_on": a.get("reviewed_on") or None}} for a in assays],
        "kits": [{
            "kit_id": k["kit_id"], "vendor": k["vendor"], "product": k["product"], "product_code": k["product_code"] or None,
            "applies_to_entries": lst(k["applies_to_entries"]), "target_scope": k["target_scope"], "format": k["format"],
            "readout": k["readout"], "url": k["url"] or None,
            "verification": {"status": k["status"], "last_checked": k.get("last_checked") or None, "source": k.get("source", "")}} for k in kits],
        "usage_stats": {"columns": USAGE_COLUMNS, "rows": stats,
                        "source_columns": ["name", "query", "retrieved_on", "method_version", "notes"],
                        "sources": sources},
        "condition_stats": conds,
        "citations": cites,
        # [label, pattern]: which words in an assay's curated readouts name each usage readout label.
        # The page compiles each pattern with new RegExp(pattern, "i").
        "readout_labels": readout_labels,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    for m in WARNINGS: print("WARNING:", m)
    for m in ERRORS: print("ERROR:  ", m)
    if ERRORS:
        print("\nNot built: fix the %d error(s) above." % len(ERRORS)); sys.exit(1)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, separators=(",", ":"))      # compact: the usage table is large
    size_mb = os.path.getsize(args.out) / 1e6
    if size_mb > MAX_JSON_MB:
        print("WARNING: %s is %.1f MB (limit we aim for: %d MB). The page will load slowly on phones; see AUDIT.md."
              % (os.path.relpath(args.out, ROOT), size_mb, MAX_JSON_MB))
    print("Built %s: %d assays, %d kits, %d usage rows, %d condition rows, %d warning(s)." % (os.path.relpath(args.out, ROOT), len(assays), len(kits), len(stats), len(conds), len(WARNINGS)))


if __name__ == "__main__":
    main()
