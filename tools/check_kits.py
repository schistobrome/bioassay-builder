#!/usr/bin/env python3
"""Check that every vendor link in data/kits.csv still opens, and print the page title to compare.

    python tools/check_kits.py
Writes data/kit_check_report.csv (kit_id, http status, page title).
Some vendor sites block scripts; a "blocked" result means needs manual review.
"""
import csv, os, re, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def check(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (assay-builder link check)"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            html = r.read(200000).decode("utf-8", "ignore")
            m = re.search(r"<title[^>]*>(.*?)</title>", html, re.S | re.I)
            return r.status, re.sub(r"\s+", " ", m.group(1)).strip() if m else ""
    except Exception as e:
        return "blocked/failed", str(e)[:80]


def main():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "data", "kits.csv"), encoding="utf-8-sig")))
    out = []
    for r in rows:
        if not r["url"]:
            out.append([r["kit_id"], "no url", ""]); continue
        s, t = check(r["url"]); out.append([r["kit_id"], s, t]); print(r["kit_id"], s, t[:70])
    with open(os.path.join(ROOT, "data", "kit_check_report.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(["kit_id", "status", "page_title"]); w.writerows(out)


if __name__ == "__main__":
    main()
