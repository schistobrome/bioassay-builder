#!/usr/bin/env python3
"""
Which equipment do the public assay texts mention?  (dispensers, liquid handlers, plate washers, readers, imagers ...)

    python scripts/equipment_scan.py                 # all summarized proteins; a few minutes; never uses the internet
    python scripts/equipment_scan.py --max-proteins 500

Writes two files in build_output/:
  equipment_mentions.csv    known equipment (a seed list of brands and generic terms): how many campaigns and
                            proteins mention each, with an example sentence
  equipment_candidates.csv  phrases like "XYZ-500 reader" or "Foo dispenser" that the seed list does NOT know,
                            so new equipment can be spotted (many will be noise; look at the top of the list)
A campaign is counted once however often its text repeats a term (ChEMBL copies of PubChem deposits count once).
"""
import argparse
import collections
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_database as bd      # noqa: E402
from build_database import m1, m2   # noqa: E402

# category -> {display name: regex}. Brand names are examples of the category, not an endorsement.
SEED = {
    "bulk reagent dispenser": {
        "Multidrop / Multidrop Combi": r"multidrop|multi-drop", "WellMate": r"wellmate", "MicroFill": r"microfill",
        "FlexDrop": r"flexdrop", "Certus Flex": r"certus", "Thermo/Titertek bulk dispenser": r"titertek",
        "bulk dispenser (generic)": r"bulk (reagent )?dispens|bottle-?top dispens|reagent dispenser"},
    "acoustic / non-contact nanoliter dispenser": {
        "Echo (Labcyte/Beckman)": r"\becho\b(?= ?(5\d\d|6\d\d|acoustic|liquid|handl|dispens|550|555))|labcyte|acoustic (liquid|droplet|dispens|transfer)",
        # not plain "mosquito": mosquito cell lines and mosquito targets are common in the assay texts
        "HP D300 / Tecan D300e": r"\bd300", "Mosquito (SPT Labtech)": r"mosquito (hts|lv|hv|x1|crystal|liquid|nanoliter|dispens)|ttp ?labtech|spt ?labtech", "BioRAPTR": r"bioraptr",
        "Dragonfly": r"dragonfly", "Nanodrop Express (Innovadyne)": r"nanodrop express|innovadyne", "Tempest": r"\btempest\b", "Nanoliter dispenser (generic)": r"nano-?lit(er|re) dispens|non-?contact dispens"},
    "pin tool / pintool": {
        "pin tool (generic)": r"pin ?tool|pin[- ]transfer|pin transfer", "V&P Scientific": r"v ?& ?p scientific", "CyBio / CyBi-Well": r"cybi-?well|cybio",
        "Kalypsys": r"kalypsys"},
    "automated liquid handler": {
        # not plain "hamilton": Hamilton syringes are hand tools
        "Biomek (Beckman)": r"biomek", "Janus (PerkinElmer)": r"\bjanus\b", "Hamilton": r"hamilton (star|microlab|nimbus|vantage|robot)|microlab (star|nimbus)|\bstarlet\b",
        "Tecan Freedom EVO / Fluent": r"freedom evo|tecan (fluent|evo)", "Agilent Bravo / Vprep": r"\bbravo\b|vprep|velocity ?11",
        "Zephyr / Sciclone": r"zephyr|sciclone", "Evolution (PerkinElmer)": r"evolution p3", "robotic platform (generic)": r"robotic|liquid[- ]handl|automated pipett|automated workstation"},
    "electronic / multichannel pipette": {"multichannel pipette": r"multi-?channel pipet","electronic pipette (Viaflo, E1-ClipTip)": r"viaflo|clip-?tip|repeater pipett|electronic pipett"},
    "microplate washer / plate handling": {
        "plate washer (named explicitly)": r"plate wash|microplate wash|\bwashers?\b",
        "wash steps (text only; may be by hand)": r"washed (3|three|4|four|5|five|times|with)|wash(ing)? (step|buffer)",
        "BioTek ELx405 / EL406 / MultiFlo": r"elx ?405|el ?406|multiflo|405 ?ts|biotek wash", "plate sealer": r"plate seal|heat seal|foil seal|pierce seal",
        "plate centrifuge": r"centrifug(ed|ation)? (the )?(plate|at)|plate centrifuge", "plate shaker": r"plate shak|orbital shak|shaken|shaker",
        # not plain "harvest": "cells were harvested" is cell culture, not a filter harvester
        "incubator (CO2)": r"co2 incubator|5% co2|humidified incubator|tissue culture incubator", "cell harvester / filter": r"harvester|filtermate|filter ?plate|glass fib(er|re) filter|unifilter"},
    "plate readers / detectors": {
        # not plain "omega" (omega-3, omega-conotoxin), "synergy" (drug synergy), "victor" or "molecular devices"
        # (that company also makes FLIPR, ImageXpress and IonWorks)
        "EnVision (PerkinElmer)": r"envision", "PHERAstar / CLARIOstar (BMG)": r"pherastar|clariostar|fluostar|polarstar|bmg labtech", "ViewLux": r"viewlux",
        "SpectraMax / FlexStation (Molecular Devices)": r"spectramax|flexstation", "Synergy / Cytation (BioTek)": r"synergy ?(h1|h4|ht|htx|2|neo\d?|mx|lx)\b|cytation|biotek",
        "Tecan Infinite / Spark / Safire": r"infinite m|tecan (spark|safire|infinite)|\bspark\b|safire", "Victor (PerkinElmer)": r"victor ?(2|3|x\d?|v|light|nivo)\b|victor multilabel", "TopCount / MicroBeta / Wallac": r"topcount|microbeta|wallac|trilux",
        "FLIPR / FDSS / kinetic imaging reader": r"\bflipr\b|fdss|hamamatsu", "Octet / Biacore / SPR system": r"octet|biacore|proteon|\bspr\b",
        "scintillation counter (generic)": r"scintillation count|liquid scintillation"},
    "imagers / high-content": {
        "ImageXpress / Opera / Operetta / CellInsight / ArrayScan": r"imagexpress|opera phenix|\bopera\b|operetta|cellinsight|arrayscan|cellomics|in cell analyzer|incell",
        "Acumen / Celigo / IncuCyte": r"acumen|celigo|incucyte", "flow cytometer (generic)": r"flow cytomet|\bfacs\b|intellicyt|\bique\b|hypercyt"},
    "automated electrophysiology / label-free": {
        "IonWorks / QPatch / Patchliner / SyncroPatch": r"ionworks|qpatch|patchliner|syncropatch|automated patch|population patch",
        "xCELLigence / impedance": r"xcelligence|rtca|impedance", "Meso Scale Discovery (ECL)": r"meso scale|\bmsd\b|sector (imager|s)|electrochemiluminesc|\becl\b"},
    "other readers seen in the data": {
        "Packard Fusion / Discovery": r"packard (fusion|discovery)|fusion plate", "Analyst HT / LJL (Molecular Devices)": r"analyst (ht|gt)|ljl analyst",
        "RUBYstar / NOVOstar / POLARstar (BMG)": r"rubystar|novostar|polarstar", "EnSpire / AlphaQuest": r"enspire|alphaquest", "Mithras / Luminoskan / Genios": r"mithras|luminoskan|genios",
        # not plain "gemini" (matches "geminin") or "M5" (the M5 muscarinic receptor); not "prominence" (an English word)
        "Gemini / M5 / i3x (Molecular Devices)": r"spectramax (gemini|m5|i3x)|\bgemini (em|xs|xps)\b|\bi3x\b", "Typhoon / phosphorimager": r"typhoon|phosphorimager|storage phosphor",
        "static light scattering / FCS": r"stargazer|fluorescence correlation|clarina", "HPLC (generic)": r"\bhplc\b|\buplc\b"},
    "mass spec / other analytical": {
        "RapidFire (Agilent)": r"rapidfire", "LC-MS (generic)": r"lc-?ms|hplc-?ms", "Caliper / LabChip mobility shift": r"caliper|labchip", "qPCR / thermal cycler": r"real-?time pcr|qpcr|lightcycler|thermal cycler|viia|quantstudio"},
}
COMPILED = [(cat, name, re.compile(rx, re.I)) for cat, d in SEED.items() for name, rx in d.items()]
CAND = re.compile(r"\b((?:[A-Z][A-Za-z0-9\-]{2,}|[A-Z]{2,}[0-9]*)(?: [A-Z0-9][A-Za-z0-9\-]*){0,2}) "
                  r"(?:plate )?(reader|dispenser|handler|washer|instrument|workstation|imager|scanner|harvester|analyzer)\b")
STOP = {"The", "This", "These", "Assay", "Detection", "Screening", "Reaction", "Reagent", "Test", "Cell", "Biochemical", "Data", "Compound", "Primary", "Protein", "Human", "Kinase"}


def snippet(text, m, width=70):
    return text[max(0, m.start() - width):m.end() + width].replace("\n", " ").strip()


def campaign_texts(records):
    """campaign key -> the texts of ALL its records joined (a ChEMBL copy joins its PubChem original).
    Before Nov 2026 only the first record of a campaign was read, so when the short ChEMBL one-liner came
    first the full PubChem protocol (where equipment is named) was never scanned."""
    originals = bd.pubchem_originals(records)
    out = collections.OrderedDict()
    for r in records:
        out.setdefault(bd.record_campaign(r, originals), []).append(r.get("text") or "")
    return {k: " || ".join(v) for k, v in out.items()}


def scan(proteins):
    """proteins: iterable of (accession, labeled records). Returns hits, candidates and counts.
    A campaign counts once however many proteins or records mention it; n_proteins counts every protein."""
    hits = collections.defaultdict(lambda: {"camps": set(), "prots": set(), "example": ""})
    cands = collections.defaultdict(lambda: {"camps": set(), "prots": set(), "example": ""})
    seen = set()
    n_prot = n_text = 0
    for acc, records in proteins:
        n_prot += 1
        for camp, text in campaign_texts(records).items():
            if camp not in seen:
                seen.add(camp)
                n_text += len(text) > 400                   # only fuller texts can name equipment
            for cat, name, rx in COMPILED:
                m = rx.search(text)
                if m:
                    h = hits[(cat, name)]
                    h["camps"].add(camp); h["prots"].add(acc)
                    h["example"] = h["example"] or snippet(text, m)
            for m in CAND.finditer(text):
                phrase = m.group(0)
                if m.group(1).split()[0] in STOP or any(rx.search(phrase) for _, _, rx in COMPILED):
                    continue
                c = cands[phrase]
                c["camps"].add(camp); c["prots"].add(acc)
                c["example"] = c["example"] or snippet(text, m)
    return hits, cands, len(seen), n_prot, n_text


def summarized_records(max_proteins=0):
    """(accession, labeled records) for every protein build_database.py summarized."""
    files = sorted(f for f in os.listdir(bd.SUMMARY) if f.endswith(".json")) if os.path.isdir(bd.SUMMARY) else []
    n = 0
    for f in files:
        s = bd.read_json(os.path.join(bd.SUMMARY, f))
        if not s:
            continue
        if max_proteins and n >= max_proteins:
            break
        records, _ = m2.build_records(s["gene"], m1.DATA, None)
        if records:
            n += 1
            yield s["acc"], records


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-proteins", type=int, default=0, help="stop after this many proteins (0 = all)")
    args = ap.parse_args(argv)
    m1.DATA = os.path.join(bd.CACHE, "m1_data")
    hits, cands, n_camp, n_prot, n_text = scan(summarized_records(args.max_proteins))
    os.makedirs(bd.OUTDIR, exist_ok=True)
    rows = [{"category": cat, "equipment": name, "n_campaigns": len(h["camps"]),
             "share_of_campaigns": round(len(h["camps"]) / max(1, n_camp), 4), "n_proteins": len(h["prots"]), "example": h["example"]}
            for (cat, name), h in hits.items()]
    rows.sort(key=lambda r: (r["category"], -r["n_campaigns"]))
    bd.write_csv(os.path.join(bd.OUTDIR, "equipment_mentions.csv"), rows,
                 ["category", "equipment", "n_campaigns", "share_of_campaigns", "n_proteins", "example"])
    crow = sorted(({"phrase": k, "n_campaigns": len(v["camps"]), "n_proteins": len(v["prots"]), "example": v["example"]}
                   for k, v in cands.items() if len(v["camps"]) >= 2), key=lambda r: -r["n_campaigns"])[:400]
    bd.write_csv(os.path.join(bd.OUTDIR, "equipment_candidates.csv"), crow, ["phrase", "n_campaigns", "n_proteins", "example"])
    print("Scanned %d campaigns on %d proteins (%d with fuller text of 400+ characters)." % (n_camp, n_prot, n_text))
    by_cat = collections.defaultdict(set)
    for (cat, name), h in hits.items():
        by_cat[cat] |= h["camps"]
    for cat in SEED:
        print("  %-45s %6d campaigns (%.1f%%)" % (cat, len(by_cat[cat]), 100.0 * len(by_cat[cat]) / max(1, n_camp)))
    print("Wrote build_output/equipment_mentions.csv (%d rows) and equipment_candidates.csv (%d rows)." % (len(rows), len(crow)))


if __name__ == "__main__":
    main()
