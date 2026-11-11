#!/usr/bin/env python3
"""
Download files from a web "Index of /..." folder (PubChem / ChEMBL / UniProt FTP-over-HTTPS), with resume.

    # 1) look first: lists files and sizes, downloads nothing
    python scripts/download_bulk.py --url https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/ --list

    # 2) download only the files you name (comma-separated; * works)
    python scripts/download_bulk.py --url https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/ ^
        --dest data/public/pubchem/Extras --only bioassays.tsv.gz,Aid2GeneidAccessionUniProt.gz,README

Safe to stop and rerun: finished files are skipped, half-finished files continue where they stopped.
Standard library only. It never deletes anything.
"""
import argparse
import fnmatch
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UA = {"User-Agent": "bioassay-builder/0.2 (personal side project)"}


def fetch(url, headers=None, method="GET", timeout=60):
    req = urllib.request.Request(url, headers=dict(UA, **(headers or {})), method=method)
    return urllib.request.urlopen(req, timeout=timeout)


def list_folder(url):
    with fetch(url) as r:
        html = r.read().decode("utf-8", errors="replace")
    files = []
    for href in re.findall(r'href="([^"]+)"', html, flags=re.I):
        if href.startswith(("?", "#", "/", "..")) or href.endswith("/") or "://" in href:
            continue
        files.append(urllib.parse.unquote(href))
    return sorted(set(files))


def size_of(url):
    try:
        with fetch(url, method="HEAD") as r:
            return int(r.headers.get("Content-Length", "-1"))
    except (urllib.error.URLError, OSError, ValueError):
        return -1


def human(n):
    if n < 0:
        return "?"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return "%.1f %s" % (n, unit)
        n /= 1024.0


def download(url, path, total):
    if total >= 0 and os.path.exists(path) and os.path.getsize(path) == total:
        print("  have  %s" % os.path.basename(path))
        return True
    part = path + ".part"
    for attempt in range(1, 8):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        try:
            with fetch(url, {"Range": "bytes=%d-" % have} if have else None, timeout=90) as r:
                mode = "ab" if (have and r.status == 206) else "wb"
                done = have if mode == "ab" else 0
                t0, last = time.time(), time.time()
                with open(part, mode) as out:
                    while True:
                        block = r.read(1 << 20)
                        if not block:
                            break
                        out.write(block)
                        done += len(block)
                        if time.time() - last > 5:
                            print("\r  %s  %s / %s" % (os.path.basename(path), human(done), human(total)), end="", flush=True)
                            last = time.time()
            if total >= 0 and os.path.getsize(part) != total:
                raise OSError("incomplete (%d of %d bytes)" % (os.path.getsize(part), total))
            os.replace(part, path)
            print("\r  done  %s (%s)%s" % (os.path.basename(path), human(os.path.getsize(path)), " " * 20))
            return True
        except urllib.error.HTTPError as e:
            if e.code == 416 and os.path.exists(part):          # already complete
                os.replace(part, path)
                return True
            err = e
        except (urllib.error.URLError, OSError) as e:
            err = e
        wait = min(60, 5 * attempt)
        print("\n  problem (%s); retrying in %d s (attempt %d of 7)" % (err, wait, attempt))
        time.sleep(wait)
    print("  GAVE UP on %s; rerun the same command later" % os.path.basename(path))
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", required=True, help="folder URL, e.g. https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/")
    ap.add_argument("--dest", default="data/public/download", help="where to save (default data/public/download)")
    ap.add_argument("--only", help="comma-separated file names or patterns, e.g. bioassays.tsv.gz,Aid2*.gz")
    ap.add_argument("--only-from", help="text file with one file name per line (made by pubchem_plan.py)")
    ap.add_argument("--skip", help="comma-separated file names or patterns to leave out")
    ap.add_argument("--list", action="store_true", help="only list files and sizes")
    ap.add_argument("--max-gb", type=float, help="refuse to download more than this many GB in total")
    args = ap.parse_args()
    url = args.url if args.url.endswith("/") else args.url + "/"
    names = list_folder(url)
    if args.only:
        pats = [p.strip() for p in args.only.split(",") if p.strip()]
        names = [n for n in names if any(fnmatch.fnmatch(n, p) for p in pats)]
    if args.only_from:
        with open(args.only_from, encoding="utf-8") as h:
            wanted = {l.strip() for l in h if l.strip()}
        missing = wanted - set(names)
        names = [n for n in names if n in wanted]
        if missing:
            print("  note: %d names in the list are not in the folder (e.g. %s)" % (len(missing), sorted(missing)[0]))
    if args.skip:
        pats = [p.strip() for p in args.skip.split(",") if p.strip()]
        names = [n for n in names if not any(fnmatch.fnmatch(n, p) for p in pats)]
    if not names:
        sys.exit("No files matched in %s" % url)
    sizes = {n: size_of(url + urllib.parse.quote(n)) for n in names}
    total = sum(s for s in sizes.values() if s > 0)
    print("%d files, %s in total" % (len(names), human(total)))
    for n in names[:60]:
        print("  %-50s %s" % (n, human(sizes[n])))
    if len(names) > 60:
        print("  ... and %d more" % (len(names) - 60))
    if args.list:
        return
    if args.max_gb and total > args.max_gb * 1024 ** 3:
        sys.exit("That is more than --max-gb %.1f. Narrow it with --only / --skip." % args.max_gb)
    os.makedirs(args.dest, exist_ok=True)
    ok = sum(download(url + urllib.parse.quote(n), os.path.join(args.dest, n), sizes[n]) for n in names)
    print("Finished: %d of %d files complete in %s" % (ok, len(names), args.dest))


if __name__ == "__main__":
    main()
