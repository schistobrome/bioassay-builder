#!/usr/bin/env python3
"""Turn the labeled public campaigns from m2_label.py into usage-share rows.

    python scripts/usage_from_m2.py EGFR --m2-dir "C:\\...\\AssayBuilder\\m2_data"

Reads <m2-dir>/<GENE>_labeled.csv (columns: campaign, in_v1_scope, function_labels, technology_labels) and
data/assays.csv (uses each assay's evidence_labels). Writes data/usage_stats.candidates.csv.
Review the file, then copy the rows you trust into data/usage_stats.csv. Nothing is added automatically.
Entry-level rows: share of classifiable campaigns that studied this function.
Readout-level rows: among those campaigns, share that used each detection technology
(only campaigns whose technology could be determined count in the denominator).
"""
import argparse, csv, os, re, collections
from datetime import date

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
split = lambda s: [x.strip() for x in re.split(r"[;|]", s or "") if x.strip()]


def campaigns_of(records):
    """Pool the labeled assays of one campaign (one paper, patent or PubChem deposit) into one row.
    Before Nov 2026 every assay row counted as a campaign, so a paper with 12 assays counted 12 times."""
    pooled = collections.OrderedDict()
    for r in records:
        c = pooled.setdefault(r.get("campaign") or r.get("id"), {"function_labels": set(), "technology_labels": set()})
        c["function_labels"].update(split(r.get("function_labels")))
        c["technology_labels"].update(split(r.get("technology_labels")))
    return [{k: "; ".join(sorted(v)) for k, v in c.items()} for c in pooled.values()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gene")
    ap.add_argument("--m2-dir", default=os.path.join(ROOT, "m2_data"))
    ap.add_argument("--assays", default=os.path.join(ROOT, "data", "assays.csv"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "usage_stats.candidates.csv"))
    a = ap.parse_args()
    path = os.path.join(a.m2_dir, "%s_labeled.csv" % a.gene)
    if not os.path.exists(path):
        raise SystemExit("Not found: %s (run m2_label.py for this gene first)" % path)
    with open(path, newline="", encoding="utf-8-sig") as f:
        records = [r for r in csv.DictReader(f) if r.get("in_v1_scope") == "yes"]
    camps = campaigns_of(records)
    total = len(camps)
    classifiable = [c for c in camps if split(c.get("function_labels"))]
    with open(a.assays, newline="", encoding="utf-8-sig") as f:
        assays = list(csv.DictReader(f))
    header = ["stat_id", "entry_id", "scope_type", "scope_id", "stage", "readout", "n_campaigns_total", "n_campaigns_classifiable",
              "n_campaigns_using", "n_literature_total", "n_literature_using", "source", "query", "retrieved_on", "method_version", "notes"]
    rows = []
    for e in assays:
        wanted = set(split(e.get("evidence_labels") or e.get("evidence_function")))
        if not wanted:
            continue
        using = [c for c in classifiable if wanted & set(split(c["function_labels"]))]
        if not using:
            continue
        base = [e["entry_id"], "gene", a.gene, "any"]
        tail = ["m2_label.py output (%s)" % os.path.basename(path), "functions: " + "; ".join(sorted(wanted)), date.today().isoformat(), "m2-labels", "auto-generated; review before use"]
        rows.append(["%s-%s-any" % (a.gene, e["entry_id"])] + base + [""] + [total, len(classifiable), len(using), "", ""] + tail)
        techs = [c for c in using if split(c.get("technology_labels"))]
        count = collections.Counter(t for c in techs for t in set(split(c["technology_labels"])))
        for tech, n in count.most_common():
            rows.append(["%s-%s-any-%s" % (a.gene, e["entry_id"], re.sub(r"\W+", "", tech)[:12])] + base + [tech] + [total, len(techs), n, "", ""] + tail)
    with open(a.out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
    print("%s: %d in-scope campaigns, %d with a function label. Wrote %d candidate rows to %s" % (a.gene, total, len(classifiable), len(rows), a.out))


if __name__ == "__main__":
    main()
