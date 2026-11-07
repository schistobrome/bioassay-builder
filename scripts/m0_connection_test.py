"""
M0 connection test for the Assay Strategy Planner.

What it does: checks that your computer can reach the public databases the app
will use (UniProt, PubChem, ChEMBL, and Europe PMC), and that their answers
look the way the app expects. It only READS public data. It changes nothing.

How to run (after installing Python from python.org):
    Windows:   python m0_connection_test.py
    Mac/Linux: python3 m0_connection_test.py

It prints a report and also saves it as m0_results.txt next to this file.
Paste the contents of m0_results.txt back into the chat.

Uses only Python's built-in libraries, so nothing else needs installing.
"""

import json
import os
import platform
import socket
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime

LOG = []


def say(text=""):
    """Print a line and remember it for the results file."""
    print(text)
    LOG.append(text)


HEADERS = {
    "User-Agent": "assay-planner-m0-test/0.1 (personal side project)",
    "Accept": "application/json",
}


def get_json(url, timeout=30):
    """Fetch a web address and return its JSON answer as Python data."""
    request = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def explain_error(error):
    """Turn a technical error into a plain-language line."""
    if isinstance(error, urllib.error.HTTPError):
        return "The server answered with error code %s (%s)." % (error.code, error.reason)
    if isinstance(error, urllib.error.URLError):
        reason = error.reason
        if isinstance(reason, ssl.SSLCertVerificationError):
            return (
                "Security certificate problem. On a Mac, open the Python folder in "
                "Applications and double-click 'Install Certificates.command', then retry. "
                "On a work network, a company firewall may be the cause."
            )
        if isinstance(reason, socket.timeout):
            return "The request timed out."
        return "Could not connect: %s" % reason
    if isinstance(error, socket.timeout):
        return "The request timed out."
    if isinstance(error, json.JSONDecodeError):
        return "The server answered, but not with data the script could read."
    if isinstance(error, (KeyError, IndexError, TypeError)):
        return "The answer arrived but was shaped differently than expected (%s: %s)." % (
            type(error).__name__,
            error,
        )
    return "%s: %s" % (type(error).__name__, error)


def short(text, limit=110):
    text = " ".join(str(text).split())
    return text if len(text) <= limit else text[: limit - 3] + "..."


# ---------------------------------------------------------------- checks


def check_uniprot(accession, expected_gene):
    """UniProt: protein length plus domain and region annotations."""
    url = (
        "https://rest.uniprot.org/uniprotkb/%s.json"
        "?fields=accession,gene_names,length,ft_domain,ft_region" % accession
    )
    data = get_json(url)
    gene = data["genes"][0]["geneName"]["value"]
    length = data.get("sequence", {}).get("length", "unknown")
    features = data.get("features", [])
    say("    Protein %s: gene %s, length %s residues, %d domain/region annotations."
        % (accession, gene, length, len(features)))
    for feature in features[:8]:
        start = feature["location"]["start"].get("value", "?")
        end = feature["location"]["end"].get("value", "?")
        say("      %-7s %5s-%-5s %s" % (feature.get("type", "?"), start, end,
                                         short(feature.get("description", ""), 60)))
    if len(features) > 8:
        say("      ... and %d more" % (len(features) - 8))
    if gene.upper() != expected_gene:
        return False, "Expected gene %s but got %s." % (expected_gene, gene)
    if not features:
        return False, "No domain or region annotations came back."
    return True, "Got %d annotations for %s." % (len(features), gene)


def check_pubchem(symbol):
    """PubChem BioAssay: assay IDs linked to a gene, then a few assay titles."""
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/genesymbol/%s/aids/JSON"
           % symbol)
    data = get_json(url)
    aids = data["IdentifierList"]["AID"]
    say("    PubChem lists %d assays linked to gene %s." % (len(aids), symbol))
    if not aids:
        return False, "No assays returned."
    time.sleep(0.5)  # PubChem asks users to stay under 5 requests per second
    first = ",".join(str(a) for a in aids[:3])
    summary_url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/%s/summary/JSON"
                   % first)
    summaries = get_json(summary_url)["AssaySummaries"]["AssaySummary"]
    for item in summaries:
        say("      AID %s | %s | %s" % (item.get("AID", "?"),
                                        short(item.get("SourceName", "?"), 30),
                                        short(item.get("Name", "?"), 70)))
    return True, "%d assays found; titles retrieved." % len(aids)


def check_chembl(symbol, expected_accession):
    """ChEMBL: find the human protein target, then list a few of its assays."""
    query = urllib.parse.quote(symbol)
    url = "https://www.ebi.ac.uk/chembl/api/data/target/search.json?q=%s&limit=20" % query
    targets = get_json(url)["targets"]
    chosen = None
    for target in targets:
        if target.get("organism") != "Homo sapiens":
            continue
        if target.get("target_type") != "SINGLE PROTEIN":
            continue
        accessions = [c.get("accession") for c in target.get("target_components", [])]
        if expected_accession in accessions:
            chosen = target
            break
    if chosen is None:
        return False, "No human single-protein target with UniProt %s found." % expected_accession
    target_id = chosen["target_chembl_id"]
    say("    Target found: %s (%s)." % (chosen.get("pref_name", "?"), target_id))
    time.sleep(0.3)
    assay_url = ("https://www.ebi.ac.uk/chembl/api/data/assay.json?target_chembl_id=%s&limit=5"
                 % target_id)
    result = get_json(assay_url)
    total = result["page_meta"]["total_count"]
    say("    ChEMBL lists %s assays for this target. First few:" % total)
    for assay in result["assays"]:
        say("      %s | type %s | %s" % (assay.get("assay_chembl_id", "?"),
                                         assay.get("assay_type", "?"),
                                         short(assay.get("description", ""), 80)))
    return True, "%s assays found." % total


def check_europepmc():
    """Europe PMC: full-text search (needed later for assays like ADP-Glo)."""
    query = urllib.parse.quote('"ADP-Glo" AND helicase')
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
           "?query=%s&format=json&pageSize=3" % query)
    data = get_json(url)
    hits = data["hitCount"]
    say("    Europe PMC found %s papers mentioning ADP-Glo and helicase. First few:" % hits)
    for paper in data["resultList"]["result"]:
        say("      %s | %s" % (paper.get("pubYear", "?"), short(paper.get("title", ""), 85)))
    return True, "%s papers found." % hits


# ---------------------------------------------------------------- runner


def run_check(name, function, *args):
    say("")
    say("[%s]" % name)
    started = time.time()
    try:
        ok, detail = function(*args)
    except Exception as error:  # we want any failure reported, not a crash
        ok, detail = False, explain_error(error)
    seconds = time.time() - started
    say("    -> %s (%.1f s): %s" % ("PASS" if ok else "FAIL", seconds, detail))
    return name, ok, detail


def main():
    say("M0 connection test for the Assay Strategy Planner")
    say("Run at: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    say("Python: %s on %s" % (sys.version.split()[0], platform.platform()))
    if sys.version_info < (3, 8):
        say("Python 3.8 or newer is needed. Please install a newer version.")
        return

    results = [
        run_check("UniProt: WRN domains", check_uniprot, "Q14191", "WRN"),
        run_check("UniProt: EGFR domains", check_uniprot, "P00533", "EGFR"),
        run_check("PubChem BioAssay: WRN assays", check_pubchem, "WRN"),
        run_check("PubChem BioAssay: EGFR assays", check_pubchem, "EGFR"),
        run_check("ChEMBL: WRN target and assays", check_chembl, "WRN", "Q14191"),
        run_check("ChEMBL: EGFR target and assays", check_chembl, "EGFR", "P00533"),
        run_check("Europe PMC: full-text search (needed later)", check_europepmc),
    ]

    say("")
    say("SUMMARY")
    for name, ok, _detail in results:
        say("    %-4s %s" % ("PASS" if ok else "FAIL", name))
    passed = sum(1 for _n, ok, _d in results if ok)
    say("    %d of %d checks passed." % (passed, len(results)))

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "m0_results.txt")
    try:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(LOG) + "\n")
        say("")
        say("Saved this report to: %s" % path)
        say("Please paste its contents back into the chat.")
    except OSError as error:
        say("")
        say("Could not save the report file (%s). Please copy the text above instead." % error)


if __name__ == "__main__":
    main()
