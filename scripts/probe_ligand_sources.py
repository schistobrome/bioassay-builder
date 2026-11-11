#!/usr/bin/env python3
"""
Probe: what do BindingDB and three neighbouring sources say about inhibitors, activators and binding sites?

    python scripts/probe_ligand_sources.py EGFR BRD4 HDAC1 GLP1R
    python scripts/probe_ligand_sources.py EGFR --sources bindingdb         # one source only
    python scripts/probe_ligand_sources.py P00533 --cutoff 10000 --max-pdb 5

Sources (all free, no key, read-only):
  bindingdb  getLigandsByUniprots (ligands with affinity type and value) and getLigandsByPDBs
             (ligands that appear in crystal structures of the protein)
  chembl     mechanism records: action_type (INHIBITOR, ACTIVATOR, AGONIST, ANTAGONIST, ALLOSTERIC ...) and binding-site comments
  pdbe       binding-site residues for the protein's crystal structures (which residues touch the ligand)
  gtopdb     IUPHAR Guide to Pharmacology: ligand "type" and "action" (Inhibitor, Activator, Agonist ...) per target

Why four? BindingDB itself records affinity (IC50, Ki, Kd, EC50), not "inhibitor / activator" labels and not
residues. An IC50 or Ki points to inhibition; an EC50 points to activation or agonism. Residue-level and allosteric
information normally come from the structure side (PDBe) and from curated mechanism records (ChEMBL, GtoPdb).
This script shows what each one actually returns for your proteins so we can decide what to build on.

It does not need to know the exact reply layout in advance: for every list of records it finds, it prints
the field names, how often each is filled, and the common values. Paste probe_report.txt back into the chat.

Files: probe_data/<source>/...  (every reply cached; reruns call nothing; --refresh to download again)
       probe_report.txt         (the summary)
Standard library only. Read-only.
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
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CACHE = os.path.join(ROOT, "probe_data")
UA = {"User-Agent": "bioassay-builder-probe/0.1 (personal side project)", "Accept": "application/json"}
ACC = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
INHIBITION_LIKE = {"ic50", "ki", "ic 50", "inhibition", "% inhibition", "kic"}
ACTIVATION_LIKE = {"ec50", "ac50", "activation", "% activation", "potency (ec50)"}
BINDING_LIKE = {"kd", "kd1", "ka", "k_d"}
LOG = []


def say(text=""):
    print(text, flush=True)
    LOG.append(text)


# ---------------------------------------------------------------- the one network door

def http_get(url, timeout=120):
    request = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def cached(source, name, url, refresh, tries=2):
    """Return (parsed, from_cache, error). Empty reply -> parsed is None. Cached on disk."""
    path = os.path.join(CACHE, source, re.sub(r"[^A-Za-z0-9._-]", "_", name) + ".json")
    if os.path.exists(path) and not refresh:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle), True, ""
    last = ""
    for attempt in range(1, tries + 1):
        try:
            text = http_get(url)
            data = json.loads(text) if text.strip() else None
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(data, handle)
            return data, False, ""
        except urllib.error.HTTPError as error:
            last = "HTTP %s %s" % (error.code, error.reason)
            if error.code not in (429, 500, 502, 503, 504):
                break
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            last = "%s: %s" % (type(error).__name__, str(error)[:100])
        except ValueError:
            last = "reply was not JSON (first characters: %r)" % text[:80]
            break
        if attempt < tries:
            time.sleep(3 * attempt)
    return None, False, last


# ---------------------------------------------------------------- generic profiling of unknown record layouts

def find_records(data):
    """Find the biggest list of dicts inside a reply, whatever the wrapper is called."""
    best = []

    def walk(node):
        nonlocal best
        if isinstance(node, list):
            dicts = [x for x in node if isinstance(x, dict)]
            if dicts and len(dicts) > len(best):
                best = dicts
            for x in node:
                walk(x)
        elif isinstance(node, dict):
            for v in node.values():
                walk(v)

    walk(data)
    return best


def profile(records, show=8):
    """Field names, fill rates and common values. Works for any layout."""
    n = len(records)
    keys = collections.OrderedDict()
    for r in records:
        for k, v in r.items():
            keys.setdefault(k, []).append(v)
    lines = []
    for k, values in keys.items():
        filled = [v for v in values if v not in (None, "", [], {})]
        flat = [json.dumps(v) if isinstance(v, (dict, list)) else str(v) for v in filled]
        distinct = collections.Counter(flat)
        text = "    %-28s filled %5.0f%%" % (k[:28], 100.0 * len(filled) / n)
        if distinct and len(distinct) <= 12:
            text += "   values: " + ", ".join("%s (%d)" % (v[:30], c) for v, c in distinct.most_common(show))
        elif distinct:
            text += "   e.g. " + "; ".join(v[:40] for v in list(distinct)[:3])
        lines.append(text)
    return lines


def tally(records, patterns):
    """Counter of the values of fields whose names match any pattern."""
    out = {}
    for r in records:
        for k, v in r.items():
            if any(re.search(p, k, re.I) for p in patterns) and not isinstance(v, (dict, list)) and v not in (None, ""):
                out.setdefault(k, collections.Counter())[str(v)] += 1
    return out


def get_any(record, *names):
    lower = {k.lower().replace(".", "_"): v for k, v in record.items()}
    for n in names:
        if n in lower and lower[n] not in (None, ""):
            return lower[n]
    return None


# ---------------------------------------------------------------- resolve protein

def resolve(query, refresh):
    q = query.strip()
    if ACC.match(q.upper()):
        acc, gene = q.upper(), ""
    else:
        url = ("https://rest.uniprot.org/uniprotkb/search?query=%s&format=json&size=1&fields=accession,gene_primary"
               % urllib.parse.quote("gene_exact:%s AND organism_id:9606 AND reviewed:true" % q))
        data, _c, err = cached("uniprot", "search_" + q, url, refresh)
        if not data or not data.get("results"):
            return None, None, err or "no human reviewed entry"
        acc, gene = data["results"][0]["primaryAccession"], q.upper()
    entry, _c, err = cached("uniprot", "entry_" + acc, "https://rest.uniprot.org/uniprotkb/%s.json" % acc, refresh)
    if entry:
        genes = entry.get("genes") or [{}]
        gene = (genes[0].get("geneName") or {}).get("value", gene)
    return acc, (gene or acc), entry


def pdb_ids(entry):
    ids = []
    for ref in (entry or {}).get("uniProtKBCrossReferences", []):
        if ref.get("database") == "PDB":
            ids.append(ref["id"])
    return ids


# ---------------------------------------------------------------- sources

def probe_bindingdb(acc, gene, entry, args):
    say("  --- BindingDB")
    url = ("https://www.bindingdb.org/rest/getLigandsByUniprots?uniprot=%s&cutoff=%d&response=application/json"
           % (acc, args.cutoff))
    data, from_cache, err = cached("bindingdb", "uniprot_%s_cutoff%d" % (acc, args.cutoff), url, args.refresh)
    if err:
        say("    getLigandsByUniprots: FAILED (%s). BindingDB is sometimes slow; try again later." % err)
    elif data is None:
        say("    getLigandsByUniprots: empty reply (no ligands at this cutoff, or the accession is not in BindingDB)")
    else:
        recs = find_records(data)
        say("    getLigandsByUniprots (cutoff %d nM)%s: %d affinity records; top-level keys: %s"
            % (args.cutoff, " [cached]" if from_cache else "", len(recs), ", ".join(list(data.keys())[:5]) if isinstance(data, dict) else type(data).__name__))
        if recs:
            for line in profile(recs):
                say(line)
            type_key = next((k for k in recs[0] if re.search(r"affinity.?type|type", k, re.I)), None)
            if type_key:
                types = collections.Counter(str(r.get(type_key, "")).strip() for r in recs)
                inh = sum(c for t, c in types.items() if t.lower() in INHIBITION_LIKE)
                act = sum(c for t, c in types.items() if t.lower() in ACTIVATION_LIKE)
                bind = sum(c for t, c in types.items() if t.lower() in BINDING_LIKE)
                say("    Reading the affinity type: inhibition-like (IC50/Ki) %d | activation-like (EC50) %d | binding only (Kd) %d | other %d"
                    % (inh, act, bind, len(recs) - inh - act - bind))
                say("    (This is a hint from the measurement type, not a curated label.)")
            ids = {str(get_any(r, "monomerid", "monomer_id", "bdb_monomerid") or "") for r in recs} - {""}
            say("    distinct ligands (monomerid): %d" % len(ids))
            for r in recs[:2]:
                say("    example: " + json.dumps(r)[:300])
    ids = pdb_ids(entry)[:args.max_pdb]
    if ids:
        url = ("https://www.bindingdb.org/rest/getLigandsByPDBs?pdb=%s&cutoff=%d&identity=100&response=application/json"
               % (",".join(ids), args.cutoff))
        data, from_cache, err = cached("bindingdb", "pdbs_%s_%s" % (acc, "_".join(ids)), url, args.refresh)
        if err:
            say("    getLigandsByPDBs (%s): FAILED (%s)" % (",".join(ids), err))
        elif data is None:
            say("    getLigandsByPDBs (%s): empty reply" % ",".join(ids))
        else:
            recs = find_records(data)
            say("    getLigandsByPDBs (%s)%s: %d records" % (",".join(ids), " [cached]" if from_cache else "", len(recs)))
            if recs:
                say("    field names: " + ", ".join(recs[0].keys()))
    else:
        say("    (UniProt lists no PDB structures, so getLigandsByPDBs was not called)")


def probe_chembl(acc, gene, entry, args):
    say("  --- ChEMBL mechanism records")
    url = ("https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession=%s&limit=5" % acc)
    data, _c, err = cached("chembl", "target_" + acc, url, args.refresh)
    if err or not data or not data.get("targets"):
        say("    no ChEMBL target found (%s)" % (err or "empty"))
        return
    targets = [t for t in data["targets"] if t.get("target_type") == "SINGLE PROTEIN"] or data["targets"]
    tid = targets[0]["target_chembl_id"]
    url = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json?target_chembl_id=%s&limit=500" % tid
    mech, from_cache, err = cached("chembl", "mechanism_" + tid, url, args.refresh)
    if err or mech is None:
        say("    mechanism: FAILED (%s)" % err)
        return
    recs = mech.get("mechanisms", [])
    say("    target %s: %d mechanism records (approved and clinical drugs only)%s" % (tid, len(recs), " [cached]" if from_cache else ""))
    if recs:
        for k, c in tally(recs, [r"^action_type$"]).items():
            say("    %s: %s" % (k, ", ".join("%s %d" % kv for kv in c.most_common())))
        for line in profile(recs, show=6)[:14]:
            say(line)
        comments = [r.get("binding_site_comment") for r in recs if r.get("binding_site_comment")]
        say("    binding_site_comment filled for %d of %d" % (len(comments), len(recs)))
        for c in comments[:3]:
            say("      e.g. " + str(c)[:160])
    else:
        say("    (no mechanism records: ChEMBL only curates mechanisms for approved and clinical compounds)")


def probe_pdbe(acc, gene, entry, args):
    say("  --- PDBe binding-site residues")
    ids = pdb_ids(entry)
    if not ids:
        say("    UniProt lists no PDB structures")
        return
    say("    %d PDB entries listed in UniProt; sampling %d" % (len(ids), min(len(ids), args.max_pdb)))
    for pdb in ids[:args.max_pdb]:
        p = pdb.lower()
        data, from_cache, err = cached("pdbe", "sites_" + p, "https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/%s" % p, args.refresh)
        if err or not data:
            say("    %s: no binding-site data (%s)" % (pdb, err or "empty"))
            continue
        sites = data.get(p) or []
        counts = []
        for s in sites:
            res = s.get("site_residues") or []
            lig = s.get("ligand_residues") or []
            counts.append((s.get("site_id", "?"), len(res), "; ".join(sorted({str(x.get("chem_comp_id", "?")) for x in lig}))[:40]))
        say("    %s: %d sites%s" % (pdb, len(sites), " [cached]" if from_cache else ""))
        for sid, n, lig in counts[:3]:
            say("      site %s: %d residues near ligand(s) %s" % (sid, n, lig))
        if sites and sites[0].get("site_residues"):
            say("      residue fields: " + ", ".join(sites[0]["site_residues"][0].keys()))
            say("      NOTE: residue numbers are the structure's own numbering; mapping them to UniProt positions needs the SIFTS mapping (next step).")
        time.sleep(0.3)


def probe_gtopdb(acc, gene, entry, args):
    say("  --- IUPHAR Guide to Pharmacology")
    url = "https://www.guidetopharmacology.org/services/targets?accession=%s" % acc
    data, from_cache, err = cached("gtopdb", "target_" + acc, url, args.refresh)
    if err:
        say("    target lookup FAILED (%s)" % err)
        return
    if not data:
        say("    no target found (not every protein is in GtoPdb: it covers drug targets)")
        return
    t = data[0] if isinstance(data, list) else data
    tid = t.get("targetId")
    say("    target %s (%s)%s" % (tid, t.get("name", "?"), " [cached]" if from_cache else ""))
    url = "https://www.guidetopharmacology.org/services/targets/%s/interactions" % tid
    inter, from_cache, err = cached("gtopdb", "interactions_%s" % tid, url, args.refresh)
    if err or not inter:
        say("    interactions: %s" % (err or "none"))
        return
    recs = find_records(inter) or (inter if isinstance(inter, list) else [])
    say("    %d ligand interactions%s" % (len(recs), " [cached]" if from_cache else ""))
    for k, c in tally(recs, [r"^type$", r"^action$", r"affinityType", r"primaryTarget"]).items():
        say("    %s: %s" % (k, ", ".join("%s %d" % kv for kv in c.most_common(8))))
    for line in profile(recs, show=5)[:12]:
        say(line)


SOURCES = collections.OrderedDict([("bindingdb", probe_bindingdb), ("chembl", probe_chembl),
                                   ("pdbe", probe_pdbe), ("gtopdb", probe_gtopdb)])


def main():
    ap = argparse.ArgumentParser(description="Probe BindingDB, ChEMBL, PDBe and GtoPdb for inhibitor/activator and binding-site information.")
    ap.add_argument("proteins", nargs="+", help="gene symbols or UniProt accessions")
    ap.add_argument("--sources", default="bindingdb,chembl,pdbe,gtopdb")
    ap.add_argument("--cutoff", type=int, default=1000, help="BindingDB affinity cutoff in nM (default 1000; larger = more ligands, slower)")
    ap.add_argument("--max-pdb", type=int, default=3, help="structures to sample per protein")
    ap.add_argument("--refresh", action="store_true", help="download again instead of using probe_data")
    args = ap.parse_args()
    wanted = [s.strip() for s in args.sources.split(",") if s.strip()]
    say("Ligand-source probe  %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for q in args.proteins:
        acc, gene, entry = resolve(q, args.refresh)
        if not acc:
            say("\n=== %s: could not resolve (%s)" % (q, entry))
            continue
        say("\n=== %s (%s)" % (gene, acc))
        for name in wanted:
            if name not in SOURCES:
                say("  unknown source '%s'" % name)
                continue
            try:
                SOURCES[name](acc, gene, entry, args)
            except Exception as error:      # one source failing must not stop the others
                say("  --- %s: crashed (%s: %s)" % (name, type(error).__name__, error))
            time.sleep(0.5)
    path = os.path.join(ROOT, "probe_report.txt")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(LOG) + "\n")
    say("\nReport saved as %s. Please paste it back into the chat." % path)


if __name__ == "__main__":
    main()
