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
USAGE_EFFECTS = {"any", "inhibition", "activation"}
ASSAY_EFFECTS = {"", "both", "inhibition", "activation"}
SCOPES = {"cell-free", "cell-based"}
ERRORS, WARNINGS = [], []


def err(msg): ERRORS.append(msg)
def warn(msg): WARNINGS.append(msg)


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

    cond_path = os.path.join(DATA, "condition_stats.csv")        # optional, written by build_database.py
    conds = []
    if os.path.exists(cond_path):
        for u in read("condition_stats.csv"):
            if u["entry_id"] not in {a["entry_id"] for a in assays}:
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

    for k in kits:
        for e in lst(k["applies_to_entries"]):
            if e not in entry_ids: err("kits.csv %s: unknown entry '%s'" % (k["kit_id"], e))
        if not k["url"] and k["status"] != "unverified": warn("kits.csv %s: no url" % k["kit_id"])

    cite_by_doi = {c["doi"].lower(): c for c in cites if c["doi"]}
    w = config["prevalence"]["weights"]
    levels = config["prevalence"]["levels"]
    min_n = config["prevalence"]["min_campaigns_for_confidence"]
    stats = []
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
        stats.append({
            "stat_id": u["stat_id"], "entry_id": u["entry_id"],
            "scope": {"type": u["scope_type"], "id": u["scope_id"]},
            "stage": u["stage"], "effect": u.get("effect") or "any", "readout": rd or None, "plate": plate,
            "campaigns": {"total": nt, "classifiable": nc, "using": nu, "share": sc},
            "literature": {"total": lt, "using": lu, "share": sl},
            "standardness": {"score": score, "level": level_for(score, levels) if score is not None else None,
                             "confidence": ("low" if (nc is not None and nc < min_n) or nc is None else "ok")},
            "source": {"name": u["source"], "query": u["query"], "retrieved_on": u["retrieved_on"] or None, "method_version": u["method_version"] or None},
            "is_example": u["method_version"] == "example",
            "notes": u["notes"]})

    out = {
        "meta": {"schema_version": "0.1", "built_on": date.today().isoformat(),
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
            "rationale": a["rationale"], "plate_formats": lst(a.get("plate_formats", "")),
            "mechanism": {"learn_more_url": a.get("learn_more_url") or None},
            "sources": [{"doi": d, **({k: cite_by_doi[d.lower()][k] for k in ("title", "year", "journal") if cite_by_doi.get(d.lower())})} for d in lst(a.get("sources_doi", ""))],
            "review": {"status": a["status"], "reviewer": a.get("reviewer") or None, "reviewed_on": a.get("reviewed_on") or None}} for a in assays],
        "kits": [{
            "kit_id": k["kit_id"], "vendor": k["vendor"], "product": k["product"], "product_code": k["product_code"] or None,
            "applies_to_entries": lst(k["applies_to_entries"]), "target_scope": k["target_scope"], "format": k["format"],
            "readout": k["readout"], "url": k["url"] or None,
            "verification": {"status": k["status"], "last_checked": k.get("last_checked") or None, "source": k.get("source", "")}} for k in kits],
        "usage_stats": stats,
        "condition_stats": conds,
        "citations": cites,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    for m in WARNINGS: print("WARNING:", m)
    for m in ERRORS: print("ERROR:  ", m)
    if ERRORS:
        print("\nNot built: fix the %d error(s) above." % len(ERRORS)); sys.exit(1)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, separators=(",", ":"))      # compact: the usage table is large
    print("Built %s: %d assays, %d kits, %d usage rows, %d condition rows, %d warning(s)." % (os.path.relpath(args.out, ROOT), len(assays), len(kits), len(stats), len(conds), len(WARNINGS)))


if __name__ == "__main__":
    main()
