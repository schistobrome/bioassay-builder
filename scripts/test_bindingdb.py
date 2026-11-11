"""
BindingDB (and companions) test for the Assay Builder.

Question being tested: can we learn, per protein, WHAT KIND of compounds exist
(inhibitors vs activators, which site they bind) from public sources?

For each protein this script tries four sources and writes what each one gives:

  1. BindingDB     ligands with affinity type (Ki / IC50 / Kd / EC50) and value, plus paper links.
                   BindingDB's own documentation lists only SMILES, affinity types and values,
                   so it does NOT say "inhibitor" or "activator" or "allosteric" outright.
                   Affinity type is a hint only (IC50/Ki usually means inhibition, EC50 suggests activation).
  2. ChEMBL        the "mechanism" table: curated drug action per target (INHIBITOR, ACTIVATOR, AGONIST,
                   ANTAGONIST, POSITIVE/NEGATIVE ALLOSTERIC MODULATOR ...) with a binding-site comment.
                   This is the closest public source to what we want.
  3. UniProt       annotated "Binding site" / "Active site" residues (with the ligand name).
  4. PDBe          residues lining the ligand in a crystal structure (first 2 PDB entries of the protein).

Nothing is changed anywhere. Raw answers are saved in bindingdb_test/ so a rerun is instant
(delete that folder to start fresh). Uses only Python's built-in libraries.

How to run (same way as the other scripts):
    Windows:   python test_bindingdb.py
    Mac/Linux: python3 test_bindingdb.py
Pick proteins (gene names or UniProt accessions):
    python test_bindingdb.py EGFR BRD4 P00533
Options:
    --cutoff 10000     BindingDB affinity cutoff in nM (default 10000)
    --skip chembl pdbe     leave out sources

Then paste bindingdb_test_report.txt back into the chat.
"""

import argparse
import collections
import hashlib
import json
import os
import re
import socket
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "bindingdb_test")
DEFAULT = ["EGFR", "WRN", "BRD4", "HDAC1", "GLP1R", "ABL1", "CDK2", "ESR1"]
HEADERS = {"User-Agent": "assay-builder-bindingdb-test/0.1 (personal side project)", "Accept": "application/json"}
PAUSE = 0.5
LOG = []


def say(text=""):
    print(text)
    LOG.append(text)


def fetch_text(url, tries=3, timeout=90):
    """Download text (cached on disk by URL). Returns (text, from_cache). Raises on failure."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, hashlib.sha1(url.encode()).hexdigest()[:16] + ".txt")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            return handle.read(), True
    last = None
    for attempt in range(1, tries + 1):
        try:
            request = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                text = response.read().decode("utf-8", "replace")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(text)
            time.sleep(PAUSE)
            return text, False
        except urllib.error.HTTPError as error:
            if error.code in (400, 404):      # a real answer: nothing found
                return "", False
            last = error
        except (urllib.error.URLError, socket.timeout, ssl.SSLError) as error:
            last = error
        time.sleep(3 * attempt)
    raise last


def explain(error):
    if isinstance(error, urllib.error.HTTPError):
        return "server error code %s (%s)" % (error.code, error.reason)
    if isinstance(error, urllib.error.URLError):
        return "could not connect (%s)" % error.reason
    return "%s: %s" % (type(error).__name__, error)


def fetch_json(url):
    text, cached = fetch_text(url)
    if not text.strip():
        return None, cached
    try:
        return json.loads(text), cached
    except ValueError:
        return {"_not_json": text[:300]}, cached


def resolve(query):
    q = query.strip()
    if re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", q.upper()):
        return q.upper()
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": "gene_exact:%s AND organism_id:9606 AND reviewed:true" % q, "format": "json", "size": 1, "fields": "accession"})
    data, _ = fetch_json(url)
    results = (data or {}).get("results") or []
    return results[0]["primaryAccession"] if results else None


# ---------------------------------------------------------------- BindingDB

def walk_records(obj, found=None):
    """Find every dict that looks like one ligand row, wherever it sits in the JSON."""
    found = [] if found is None else found
    if isinstance(obj, dict):
        keys = " ".join(obj.keys()).lower()
        if "affinity" in keys and ("smile" in keys or "monomer" in keys):
            found.append(obj)
        else:
            for value in obj.values():
                walk_records(value, found)
    elif isinstance(obj, list):
        for item in obj:
            walk_records(item, found)
    return found


def pick(row, *needles):
    for key, value in row.items():
        low = key.lower()
        if all(n in low for n in needles):
            return value
    return None


def value_of(row):
    for key, value in row.items():
        if key.lower() == "affinity":
            return value
    for key, value in row.items():
        if "affinity" in key.lower() and "type" not in key.lower():
            return value
    return None


def to_number(value):
    try:
        return float(re.sub(r"[<>~= ]", "", str(value)))
    except (TypeError, ValueError):
        return None


def do_bindingdb(acc, cutoff):
    say("  BindingDB:")
    urls = [
        "https://bindingdb.org/rest/getLigandsByUniprots?uniprot=%s&cutoff=%d&response=application/json" % (acc, cutoff),
        "https://bindingdb.org/rest/getLigandsByUniprot?uniprot=%s;%d&response=application/json" % (acc, cutoff),
    ]
    rows, used = [], None
    for url in urls:
        try:
            data, cached = fetch_json(url)
        except Exception as error:
            say("    %s -> FAILED, %s" % (url.split("?")[0].rsplit("/", 1)[1], explain(error)))
            continue
        if not data:
            say("    %s -> empty answer (BindingDB says no match)" % url.split("?")[0].rsplit("/", 1)[1])
            continue
        if "_not_json" in data:
            say("    %s -> answer is not JSON, starts: %s" % (url.split("?")[0].rsplit("/", 1)[1], data["_not_json"][:120]))
            continue
        rows = walk_records(data)
        used = url
        say("    %s -> %d ligand rows%s" % (url.split("?")[0].rsplit("/", 1)[1], len(rows), " (cached)" if cached else ""))
        if rows:
            break
    if not rows:
        return
    say("    Field names in a row: " + ", ".join(sorted(rows[0].keys())))
    types = collections.Counter(str(pick(r, "affinity", "type") or pick(r, "type") or "?") for r in rows)
    say("    Affinity types: " + ", ".join("%s %d" % kv for kv in types.most_common()))
    mono = {str(pick(r, "monomer") or pick(r, "smile")) for r in rows}
    say("    Distinct compounds: %d" % len(mono))
    values = [to_number(value_of(r)) for r in rows]
    values = [v for v in values if v is not None]
    if values:
        say("    Affinity values (nM): best %.3g, median %.3g; %d rows at or below 100 nM"
            % (min(v for v in values if v > 0), sorted(values)[len(values) // 2], sum(1 for v in values if v <= 100)))
    papers = {str(pick(r, "doi") or pick(r, "pmid")) for r in rows if pick(r, "doi") or pick(r, "pmid")}
    say("    Distinct papers (DOI/PMID) referenced: %d" % len(papers))
    ec = types.get("EC50", 0) + types.get("Potency", 0)
    inh = types.get("IC50", 0) + types.get("Ki", 0)
    say("    Hint only: %d IC50/Ki rows (inhibition-type) vs %d EC50/Potency rows (activation-type). "
        "BindingDB does not state inhibitor/activator or binding site." % (inh, ec))


# ------------------------------------------------------------------ ChEMBL

def do_chembl(acc):
    say("  ChEMBL mechanisms:")
    base = "https://www.ebi.ac.uk/chembl/api/data/"
    try:
        data, _ = fetch_json(base + "target.json?target_components__accession=%s&target_type=SINGLE%%20PROTEIN&limit=5" % acc)
    except Exception as error:
        say("    target lookup FAILED, %s" % explain(error))
        return
    targets = (data or {}).get("targets") or []
    if not targets:
        say("    no ChEMBL single-protein target for %s" % acc)
        return
    tid = targets[0]["target_chembl_id"]
    say("    target %s (%s)" % (tid, targets[0].get("pref_name", "")))
    try:
        mech, _ = fetch_json(base + "mechanism.json?target_chembl_id=%s&limit=1000" % tid)
    except Exception as error:
        say("    mechanism lookup FAILED, %s" % explain(error))
        return
    items = (mech or {}).get("mechanisms") or []
    say("    %d mechanism rows" % len(items))
    if not items:
        say("    (ChEMBL has no curated drug mechanism for this target; it may have bioactivity data only)")
        return
    say("    Action types: " + ", ".join("%s %d" % kv for kv in collections.Counter(
        i.get("action_type") or "?" for i in items).most_common()))
    sites = [i for i in items if i.get("binding_site_comment")]
    say("    Rows with a binding-site comment: %d" % len(sites))
    for i in sites[:3]:
        say("      %s | %s | %s" % (i.get("action_type"), i.get("mechanism_of_action"), str(i.get("binding_site_comment"))[:110]))
    say("    Row fields: " + ", ".join(sorted(items[0].keys())))


# ----------------------------------------------------------------- UniProt

def do_uniprot(acc):
    say("  UniProt binding and active sites:")
    try:
        entry, _ = fetch_json("https://rest.uniprot.org/uniprotkb/%s.json" % acc)
    except Exception as error:
        say("    FAILED, %s" % explain(error))
        return []
    feats = [f for f in (entry or {}).get("features", []) if f.get("type") in ("Binding site", "Active site", "Site")]
    count = collections.Counter(f["type"] for f in feats)
    say("    " + (", ".join("%s %d" % kv for kv in count.most_common()) if feats else "none annotated"))
    ligands = collections.Counter((f.get("ligand") or {}).get("name", "") for f in feats if f.get("ligand"))
    if ligands:
        say("    Ligands named: " + ", ".join("%s (%d)" % kv for kv in ligands.most_common(6)))
    for f in feats[:3]:
        loc = f.get("location", {})
        say("      %s %s-%s %s" % (f["type"], (loc.get("start") or {}).get("value"), (loc.get("end") or {}).get("value"),
                                 (f.get("ligand") or {}).get("name", f.get("description", ""))))
    pdbs = [r["id"] for r in (entry or {}).get("uniProtKBCrossReferences", []) if r.get("database") == "PDB"]
    say("    PDB entries listed: %d" % len(pdbs))
    return pdbs


# -------------------------------------------------------------------- PDBe

def do_pdbe(pdbs):
    say("  PDBe ligand-binding residues (first 2 PDB entries):")
    for pdb in pdbs[:2]:
        try:
            data, _ = fetch_json("https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/%s" % pdb.lower())
        except Exception as error:
            say("    %s FAILED, %s" % (pdb, explain(error)))
            continue
        sites = (data or {}).get(pdb.lower()) or []
        say("    %s: %d binding sites%s" % (pdb, len(sites), "" if sites else " (nothing returned)"))
        for s in sites[:2]:
            residues = s.get("site_residues") or []
            say("      site %s | %s | %d residues, e.g. %s" % (
                s.get("site_id", "?"), str(s.get("details", ""))[:60], len(residues),
                ", ".join("%s%s" % (r.get("residue_number", ""), r.get("chem_comp_id", "")) for r in residues[:6])))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("proteins", nargs="*", help="gene names or UniProt accessions")
    parser.add_argument("--cutoff", type=int, default=10000)
    parser.add_argument("--skip", nargs="*", default=[], choices=["bindingdb", "chembl", "uniprot", "pdbe"])
    args = parser.parse_args()
    proteins = args.proteins or DEFAULT

    say("BindingDB and companion-source test for the Assay Builder")
    say("Run at: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    say("Affinity cutoff: %d nM. Raw answers are saved in %s" % (args.cutoff, CACHE))
    for query in proteins:
        say("")
        say("=== %s ===" % query)
        try:
            acc = resolve(query)
        except Exception as error:
            say("  UniProt lookup FAILED, %s" % explain(error))
            continue
        if not acc:
            say("  no human reviewed UniProt entry found")
            continue
        say("  UniProt accession: %s" % acc)
        if "bindingdb" not in args.skip:
            try:
                do_bindingdb(acc, args.cutoff)
            except Exception as error:
                say("  BindingDB step crashed: %s: %s" % (type(error).__name__, error))
        if "chembl" not in args.skip:
            try:
                do_chembl(acc)
            except Exception as error:
                say("  ChEMBL step crashed: %s: %s" % (type(error).__name__, error))
        pdbs = []
        if "uniprot" not in args.skip:
            try:
                pdbs = do_uniprot(acc)
            except Exception as error:
                say("  UniProt step crashed: %s: %s" % (type(error).__name__, error))
        if "pdbe" not in args.skip and pdbs:
            try:
                do_pdbe(pdbs)
            except Exception as error:
                say("  PDBe step crashed: %s: %s" % (type(error).__name__, error))
    say("")
    say("DONE.")
    path = os.path.join(HERE, "bindingdb_test_report.txt")
    try:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(LOG) + "\n")
        say("Report saved as: %s (please paste it back into the chat)" % path)
    except OSError as error:
        say("Could not save the report (%s). Please copy the text above." % error)


if __name__ == "__main__":
    main()
