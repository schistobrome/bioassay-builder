#!/usr/bin/env python
"""Find assay types and activity directions that the literature shows for a protein but our campaign
data does not.

This answers the coverage question, not the methodology question: if HTRF is used for protein A in
papers, but no ChEMBL or PubChem campaign for protein A states HTRF, the app will never list HTRF for
protein A. That is a gap in what we present, and it comes from how the database was assembled rather
than from biology.

Three gaps are reported per protein:

  readout gap    Technology labels found in paper full text that no campaign for that protein states.
                 Detection uses the SAME term table (scripts/term_table.csv) and the same normalisation
                 as scripts/m2_label.py, so a readout named here is one the pipeline would have
                 labelled had the text been in a campaign description.

  direction      NOT reported. Keyword detection of activity direction was tried and does not work:
                 every activation claim it produced was a false positive, because "agonist" usually
                 belongs to another protein in the same sentence, "potentiator" means a drug
                 combination, "activator" is often part of a protein name, and "PAM" means protospacer
                 adjacent motif in any paper with a CRISPR section. Direction must come from the
                 structured extraction (effect_direction and pharmacology in
                 schema/paper_extraction.schema.json), not from keywords.

  entry gap      Playbook entries that would gain protein-specific evidence if the paper readouts were
                 counted, i.e. entries the app currently shows for this protein on class evidence only.

Inputs are already on disk: build_cache/papers/search/<GENE>.json for which papers belong to which
protein, and build_cache/papers/fulltext/<PMCID>.xml for the text. Nothing is downloaded.

Examples
--------
  py scripts\\readout_gaps.py
  py scripts\\readout_gaps.py --genes CASP4,GPBAR1 --verbose
  py scripts\\readout_gaps.py --min-papers 3

Output
------
  build_output/papers/readout_gaps.csv   one row per (protein, readout) gap
  build_output/papers/readout_gaps.md    summary by protein and by readout
"""

import argparse
import collections
import csv
import glob
import json
import os
import re
import sys
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
CACHE = os.path.join(ROOT, "build_cache")
SEARCH_DIR = os.path.join(CACHE, "papers", "search")
FULL_DIR = os.path.join(CACHE, "papers", "fulltext")
OUT_DIR = os.path.join(ROOT, "build_output", "papers")

sys.path.insert(0, HERE)
import target_subject as tsub


def log(msg):
    print(msg)
    sys.stdout.flush()


def read_csv_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_json(path, default=None):
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return default


def normalise(text):
    """Same idea as m2_label.normT: lower case, punctuation to spaces, so 'TR-FRET' and 'tr fret' match."""
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", (text or "").lower()))


def load_technology_terms():
    """(term, compiled pattern) for every technology term in the term table, built from its own
    synonyms column so detection matches what the pipeline would have labelled."""
    out = []
    for row in read_csv_rows(os.path.join(HERE, "term_table.csv")):
        if (row.get("axis") or "").strip().lower() != "technology":
            continue
        term = (row.get("term") or "").strip()
        if not term:
            continue
        alts = [normalise(term)]
        for syn in (row.get("synonyms") or "").split("|"):
            syn = normalise(syn)
            if syn:
                alts.append(syn)
        alts = sorted({a for a in alts if len(a) >= 3}, key=len, reverse=True)
        if alts:
            # Word boundaries are essential. Without them "spr" matches inside "sprayed" and "cetsa"
            # inside "acetsalicylate", which made SPR look universal in a first run.
            out.append((term, re.compile(r"\b(?:%s)\b" % "|".join(re.escape(a) for a in alts))))
    return out


# Techniques that appear in the methods of almost any cell-biology paper, whether or not they are an
# assay for the target. Reported separately so they cannot dominate the gap list.
# Synonyms in the term table that are also ordinary laboratory words. "caliper" is listed as a
# synonym of the Caliper / LabChip mobility-shift platform, and in a methods section it far more often
# means the vernier caliper used to measure tumour volume: "tumor dimensions were assessed twice
# weekly with a digital caliper". Three of the MTAP false positives were exactly that sentence.
# Guarded here rather than fixed in scripts/term_table.csv, because editing the term table forces a
# full m2 relabel. If the owner wants it fixed at source, the synonym should become "caliper assay|
# caliper labchip" and this guard can go.
AMBIGUOUS_SYNONYM = {
    "Mobility shift (Caliper / LabChip)": re.compile(
        r"mobility shift|labchip|caliper (?:assay|instrument|platform|microfluidic)|"
        r"\b(?:lc3000|ez reader)\b"),
}


def readout_confirmed(term, norm_text):
    """For a term with an ambiguous synonym, require the unambiguous wording somewhere in the text."""
    guard = AMBIGUOUS_SYNONYM.get(term)
    return True if guard is None else bool(guard.search(norm_text))


GENERIC_TECHNIQUES = {
    "Western blot", "Gel-based (PAGE / EMSA)", "Absorbance / colorimetric",
    "High-content imaging", "Mass spectrometry", "Crystallography / cryo-EM",
}


# Activation and inhibition language. Deliberately narrow, because these decide a direction claim.
# EC50 and Emax are NOT activation markers: both are routinely reported for inhibitors in cell assays,
# and including EC50 made inhibition-only proteins such as WRN appear to have activators.
# "activation" alone is not a marker either: "caspase-4 activation" and "inflammasome activation"
# describe the biology being measured, not a compound that activates the target.
# Tightened after validation. The loose version produced five activation claims for WRN, all false:
#   "E1 ubiquitin-activating enzyme"   - a protein name
#   "activator protein 1 (AP-1)"       - a protein name
#   "potentiated with irinotecan"      - drug combination, not target activation
#   "agonistic factor for tumorigenesis" - general biology
# So: nouns only, no "-ing" or "-ic" forms, and the phrase must also pass the proximity test.
ACTIVATION_RE = re.compile(
    r"\bagonists?\b|\bagonism\b|\bpotentiators?\b|positive allosteric modulators?\b|\bPAMs?\b|"
    r"\bactivators? of\b|\ballosteric activators?\b")
ACTIVATION_EXCLUDE_RE = re.compile(
    r"activating enzyme|activator protein|\bAP-1\b|plasminogen activator|"
    r"activator of transcription|co-?activator")
INHIBITION_RE = re.compile(
    r"\binhibitors?\b|\binhibition\b|\bantagonists?\b|\bantagonism\b|\bIC50\b|\bblockers?\b|"
    r"\bdegraders?\b|negative allosteric modulators?\b|\bNAMs?\b")


def direction_signal(whole_text, alias_pattern, window=600):
    """(activation, inhibition) booleans for one paper. Both require the wording to sit near a mention
    of the target, so a reference title or an unrelated pathway cannot set the direction."""
    norm_alias = alias_pattern
    def near(pattern, exclude=None):
        for m in pattern.finditer(whole_text):
            seg = whole_text[max(0, m.start() - 90):m.start() + 90]
            if exclude is not None and exclude.search(seg):
                continue
            if norm_alias is None:
                return True
            lo = max(0, m.start() - window)
            if norm_alias.search(normalise(whole_text[lo:m.start() + window])):
                return True
        return False
    return near(ACTIVATION_RE, ACTIVATION_EXCLUDE_RE), near(INHIBITION_RE)


def target_aliases(gene, acc):
    """Normalised strings that mean "this protein", for the proximity test."""
    alts = {normalise(gene)}
    entry = read_json(os.path.join(CACHE, "m1_data", "uniprot_%s.json" % acc), {}) or {}
    desc = entry.get("proteinDescription") or {}
    rec = ((desc.get("recommendedName") or {}).get("fullName") or {}).get("value")
    if rec:
        alts.add(normalise(rec))
    for alt in (desc.get("alternativeNames") or [])[:3]:
        v = (alt.get("fullName") or {}).get("value")
        if v:
            alts.add(normalise(v))
    for g in entry.get("genes") or []:
        for syn in (g.get("synonyms") or [])[:4]:
            v = normalise((syn.get("value") or ""))
            if len(v) >= 3:
                alts.add(v)
    return sorted({a for a in alts if len(a) >= 3}, key=len, reverse=True)


def near_target(norm_text, readout_pattern, alias_pattern, window=600, spans=None):
    """True if the readout is named within `window` characters of a mention of the target.
    This is the difference between "the paper mentions HTRF" and "HTRF was used for this protein".

    `spans` are precomputed mention offsets, normally from target_subject.subject_positions, which
    excludes biomarker mentions: being near "MTAP-deleted HCT116 cells" is not evidence that the
    readout measures MTAP. Recomputing them per readout term would be 32x the work for no gain."""
    if alias_pattern is None and spans is None:
        return bool(readout_pattern.search(norm_text))
    if spans is None:
        spans = [m.start() for m in alias_pattern.finditer(norm_text)]
    if not spans:
        return False
    for m in readout_pattern.finditer(norm_text):
        pos = m.start()
        for s in spans:
            if abs(s - pos) <= window:
                return True
    return False

METHODS_RE = re.compile(r"method|experimental|materials|protocol", re.I)


def element_text(el):
    return " ".join("".join(el.itertext()).split())


def paper_text(pmcid, methods_only=True):
    """Full text, or just the methods sections. Methods-only is the stricter test: a readout named in
    a methods section was actually used, while one named in an introduction may only be cited."""
    path = os.path.join(FULL_DIR, pmcid + ".xml")
    if not os.path.exists(path):
        return "", ""
    try:
        with open(path, "rb") as handle:
            root = ET.fromstring(handle.read())
    except (ET.ParseError, OSError):
        return "", ""
    body = root.find(".//body")
    if body is None:
        return "", ""
    whole = element_text(body)
    if not methods_only:
        return whole, whole
    chunks = []
    for sec in body.iter("sec"):
        title_el = sec.find("title")
        title = element_text(title_el) if title_el is not None else ""
        if METHODS_RE.search(title):
            chunks.append(element_text(sec))
    return whole, (" ".join(chunks) if chunks else whole)


def campaign_data():
    """Per gene: readouts and directions already in our campaign data, plus the classes and counts."""
    out = {}
    for path in glob.glob(os.path.join(CACHE, "summary", "*.json")):
        s = read_json(path)
        if not s:
            continue
        gene = (s.get("scope") or s.get("gene") or "").strip()
        if not gene:
            continue
        camps = s.get("campaigns") or []
        classifiable = [c for c in camps if c.get("f")]
        readouts = collections.Counter()
        directions = collections.Counter()
        for c in classifiable:
            for t in c.get("t") or []:
                readouts[t] += 1
            directions[c.get("e") or "unspecified"] += 1
        out[gene] = {"acc": s.get("acc"), "classes": s.get("classes") or [],
                     "n_campaigns": len(camps), "n_classifiable": len(classifiable),
                     "readouts": readouts, "directions": directions,
                     "n_stating_readout": sum(1 for c in classifiable if c.get("t"))}
    return out


def playbook():
    """Entries with their classes, evidence labels and effect restriction."""
    rows = []
    for a in read_csv_rows(os.path.join(DATA, "assays.csv")):
        classes = {x.strip() for x in (a.get("classes") or "").split(";") if x.strip()}
        labels = set()
        for part in (a.get("evidence_labels") or a.get("evidence_function") or "").split(";"):
            for x in part.split("|"):
                if x.strip():
                    labels.add(x.strip())
        readouts = [x.strip() for x in (a.get("readouts") or "").split(";") if x.strip()]
        rows.append({"entry_id": a["entry_id"], "name": a.get("assay_name") or "",
                     "classes": classes, "labels": labels, "readouts": readouts,
                     "effect": (a.get("effect") or "").strip()})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", help="comma-separated subset; default every gene with cached papers")
    ap.add_argument("--min-papers", type=int, default=2,
                    help="a readout must appear in at least this many of the protein's papers (default 2); "
                         "one paper is an anecdote")
    ap.add_argument("--window", type=int, default=600,
                    help="a readout counts only if named within this many characters of a mention of "
                         "the target (default 600, roughly a paragraph); 0 disables the proximity test")
    ap.add_argument("--whole-text", action="store_true",
                    help="search the whole body rather than the methods sections only (looser, noisier)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    terms = load_technology_terms()
    if not terms:
        sys.exit("No technology terms: scripts/term_table.csv not readable.")
    camp = campaign_data()
    entries = playbook()
    want = {g.strip().upper() for g in args.genes.split(",")} if args.genes else None

    os.makedirs(OUT_DIR, exist_ok=True)
    gap_rows = []
    per_gene = {}

    for spath in sorted(glob.glob(os.path.join(SEARCH_DIR, "*.json"))):
        d = read_json(spath) or {}
        gene = d.get("gene")
        if not gene or (want and gene.upper() not in want):
            continue
        pmcids = [r.get("pmcid") for r in (d.get("results") or []) if r.get("pmcid")]
        pmcids = [p for p in pmcids if os.path.exists(os.path.join(FULL_DIR, p + ".xml"))]
        if not pmcids:
            continue
        base = camp.get(gene)
        if not base:
            continue

        # Family-name aliases are dropped before anything else: UniProt lists "Purinergic receptor"
        # under P2RX3 and "CK1" under CSNK1A1L, and either one makes every sibling's paper look like
        # this protein's. The gene symbol is always kept, however short.
        aliases = tsub.specific_aliases(target_aliases(gene, base["acc"]), normalise(gene))
        alias_pattern = re.compile(r"\b(?:%s)\b" % "|".join(re.escape(a) for a in aliases)) if aliases else None
        dropped = sorted(set(target_aliases(gene, base["acc"])) - set(aliases))
        if dropped:
            log("  %-10s dropped family alias(es): %s" % (gene, ", ".join(dropped)))

        found = collections.Counter()          # readout -> papers naming it NEAR the target, in methods
        act_papers = inh_papers = 0
        n_read = 0
        for pmcid in pmcids:
            whole, methods = paper_text(pmcid, methods_only=not args.whole_text)
            if not methods:
                continue
            n_read += 1
            norm = normalise(methods)
            # Subject mentions only. A readout sitting beside "MTAP-null cells" is measuring
            # whatever the paper's real target is, not MTAP.
            spans = tsub.subject_positions(norm, alias_pattern)
            seen = set()
            for term, pattern in terms:
                if near_target(norm, pattern, alias_pattern, args.window, spans) \
                        and readout_confirmed(term, norm):
                    seen.add(term)
            for term in seen:
                found[term] += 1
            is_act, is_inh = direction_signal(whole, alias_pattern, args.window)
            act_papers += 1 if is_act else 0
            inh_papers += 1 if is_inh else 0

        have = set(base["readouts"])
        new_all = {t: n for t, n in found.items() if t not in have and n >= args.min_papers}
        new = {t: n for t, n in new_all.items() if t not in GENERIC_TECHNIQUES}
        new_generic = {t: n for t, n in new_all.items() if t in GENERIC_TECHNIQUES}
        confirmed = {t: n for t, n in found.items() if t in have}

        # directions our campaigns record
        dirs = base["directions"]
        camp_act = dirs.get("activation", 0) + dirs.get("mixed", 0)
        camp_inh = dirs.get("inhibition", 0) + dirs.get("mixed", 0)
        # Direction is computed but NOT reported as a gap. Validation found every keyword-derived
        # activation claim to be a false positive, across WRN, CASP4 and CSNK1A1L:
        #   "TLR4 agonist", "NLRP3 agonist"         - an agonist of a different protein nearby
        #   "potentiators of existing antifungals"  - drug combination
        #   "agonist-challenged platelets"          - a platelet agonist
        #   "PAM is highlighted in bold"            - protospacer adjacent motif, in a CRISPR section
        #   "activator of Casein kinase 1 epsilon"  - a different family member
        # The vocabulary is too overloaded for keyword matching to settle direction. Direction has to
        # come from the structured extraction, where a model reads the sentence and fills
        # effect_direction and pharmacology per assay. Kept here only as a diagnostic counter.
        direction_gap = ""

        # which entries would gain protein-specific evidence from the new readouts
        cls = set(base["classes"]) | {"generic"}
        gained = []
        for e in entries:
            if not (e["classes"] & cls):
                continue
            if e["labels"] & set():      # placeholder, evidence labels are function not technology
                pass
            for rd in e["readouts"]:
                nrd = normalise(rd)
                for t in new:
                    if normalise(t) in nrd or nrd in normalise(t):
                        gained.append(e["entry_id"])
                        break
                else:
                    continue
                break
        gained = sorted(set(gained))

        per_gene[gene] = {"papers": n_read, "classes": ";".join(base["classes"]),
                          "n_classifiable": base["n_classifiable"],
                          "n_stating_readout": base["n_stating_readout"],
                          "have": sorted(have), "new": new, "new_generic": new_generic,
                          "confirmed": confirmed,
                          "direction_gap": direction_gap, "act_papers": act_papers,
                          "inh_papers": inh_papers, "camp_act": camp_act, "camp_inh": camp_inh,
                          "entries_gained": gained}
        for term, n in sorted(new.items(), key=lambda kv: -kv[1]):
            gap_rows.append({"gene": gene, "classes": ";".join(base["classes"]),
                             "readout_new_in_papers": term, "papers_naming_it": n,
                             "papers_searched": n_read,
                             "campaigns_classifiable": base["n_classifiable"],
                             "campaigns_stating_any_readout": base["n_stating_readout"],
                             "readouts_already_known": "; ".join(sorted(have)),
                             "entries_that_would_gain": "; ".join(gained),
                             "direction_gap": direction_gap})

    with open(os.path.join(OUT_DIR, "readout_gaps.csv"), "w", newline="", encoding="utf-8") as handle:
        cols = ["gene", "classes", "readout_new_in_papers", "papers_naming_it", "papers_searched",
                "campaigns_classifiable", "campaigns_stating_any_readout", "readouts_already_known",
                "entries_that_would_gain", "direction_gap"]
        w = csv.DictWriter(handle, fieldnames=cols)
        w.writeheader()
        w.writerows(gap_rows)

    write_markdown(per_gene, args)
    log("proteins analysed: %d" % len(per_gene))
    log("readout gaps found: %d across %d protein(s)"
        % (len(gap_rows), len({r["gene"] for r in gap_rows})))
    log("  build_output/papers/readout_gaps.csv")
    log("  build_output/papers/readout_gaps.md")

    if args.verbose:
        for gene, g in sorted(per_gene.items()):
            log("\n=== %s (%s), %d papers read" % (gene, g["classes"] or "-", g["papers"]))
            log("  campaigns: %d classifiable, %d state a readout"
                % (g["n_classifiable"], g["n_stating_readout"]))
            log("  already known: %s" % (", ".join(g["have"]) or "NONE"))
            log("  NEW in papers: %s" % (", ".join("%s (%d papers)" % kv for kv in
                                                   sorted(g["new"].items(), key=lambda kv: -kv[1])) or "none"))
            log("  (generic techniques, reported separately: %s)"
                % (", ".join("%s (%d)" % kv for kv in
                             sorted(g["new_generic"].items(), key=lambda kv: -kv[1])) or "none"))
            log("  confirmed:     %s" % (", ".join("%s (%d)" % kv for kv in
                                                   sorted(g["confirmed"].items(), key=lambda kv: -kv[1])) or "none"))
            if g["direction_gap"]:
                log("  DIRECTION GAP: %s (papers: %d activation, %d inhibition; campaigns: %d activation, %d inhibition)"
                    % (g["direction_gap"], g["act_papers"], g["inh_papers"], g["camp_act"], g["camp_inh"]))
            if g["entries_gained"]:
                log("  entries that would gain evidence: %s" % ", ".join(g["entries_gained"]))
    return 0


def write_markdown(per_gene, args):
    lines = ["# Assay types the literature shows but our campaign data does not\n"]
    lines.append("Readouts detected in paper methods sections using the same term table "
                 "(`scripts/term_table.csv`) the pipeline labels campaigns with, then compared with the "
                 "readouts our campaigns state for the same protein. A readout is reported only if at "
                 "least %d of the protein's papers name it.\n" % args.min_papers)
    lines.append("\n## Per protein\n")
    lines.append("| protein | classes | papers | campaigns stating a readout | readouts we have | "
                 "readouts only in papers | direction gap |")
    lines.append("|---|---|---|---|---|---|---|")
    for gene, g in sorted(per_gene.items(), key=lambda kv: -len(kv[1]["new"])):
        lines.append("| **%s** | %s | %d | %d of %d | %d | **%s** | %s |" % (
            gene, (g["classes"] or "-").replace(";", ", "), g["papers"],
            g["n_stating_readout"], g["n_classifiable"], len(g["have"]),
            ", ".join(sorted(g["new"], key=lambda t: -g["new"][t])) or "none",
            g["direction_gap"] or "-"))

    counts = collections.Counter()
    for g in per_gene.values():
        for term in g["new"]:
            counts[term] += 1
    lines.append("\n## Which readouts are most often missing\n")
    lines.append("| readout | proteins where it appears only in papers |")
    lines.append("|---|---|")
    for term, n in counts.most_common():
        lines.append("| %s | %d |" % (term, n))

    gained = collections.Counter()
    for g in per_gene.values():
        for e in g["entries_gained"]:
            gained[e] += 1
    if gained:
        lines.append("\n## Playbook entries that would gain protein-specific evidence\n")
        lines.append("| entry | proteins |")
        lines.append("|---|---|")
        for e, n in gained.most_common():
            lines.append("| %s | %d |" % (e, n))

    lines.append("\n## How to read this\n")
    lines.append("- A readout in the \"only in papers\" column means the app currently cannot list it "
                 "for that protein, because no campaign for that protein states it. The readout order "
                 "on the card, and the \"most used for this protein\" line, both come from campaign "
                 "counts alone.\n")
    lines.append("- Activity direction is deliberately absent from this report. Keyword detection of "
                 "activation was tried and every claim it produced was wrong, because the vocabulary "
                 "is overloaded: \"agonist\" normally refers to another protein in the same sentence, "
                 "\"potentiator\" means a drug combination, \"activator\" is often part of a protein "
                 "name, and \"PAM\" means protospacer adjacent motif wherever a paper has a CRISPR "
                 "section. Direction has to come from the structured extraction instead.\n")
    lines.append("- Detection is keyword based on methods sections, so it tells you a readout was "
                 "named in a methods section, not that it was the paper's primary assay. Treat the "
                 "output as a list of candidates for review, not as counts to publish.\n")
    with open(os.path.join(OUT_DIR, "readout_gaps.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
