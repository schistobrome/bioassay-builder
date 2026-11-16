# How the extracted paper JSONs get into the app and the database

Companion to `PAPER_EXTRACTION.md`, which covers getting the JSON out of papers. This covers what
happens next: where each field lands, what has to be built, and the three things that can go wrong.

## The governing principle: papers are a third evidence tier, not more campaigns

The app already has two tiers behind every number: the protein's own campaigns, and its function class
when the protein has too few. Paper extractions become a **third, separate tier**. They are never
folded into the campaign counts.

That is not fastidiousness, it is correctness. The card says "measures a function used in 87% of 1,847
campaigns", and `usage_stats.csv` defines that share as *campaigns using / campaigns classifiable*.
Adding papers to the numerator would change what the sentence means while leaving the sentence alone.
Worse, we cannot currently tell whether a mined paper **is already** one of those campaigns (see
"Deduplication" below), so folding them in would double count by an unknown amount.

So: papers get their own tables, their own aggregation, and their own labelled row in the UI. Where
they feed curated columns in `data/assays.csv`, they do it the way everything else does — as a review
chunk the owner signs off.

## Where each part of the JSON goes

| JSON path | Destination | How it gets there |
|---|---|---|
| `conditions.*` (ph, dmso_pct, atp_um, volume_ul, buffer, detergent …) | `data/paper_conditions.csv` **(new)** → `docs/assay_db.json` | new aggregator, new UI tier |
| `quality.z_prime`, `signal_to_background`, `dmso_tolerance_*` | `data/assays.csv` → `qc_metric`, `dmso_pct`, plus a new `z_prime_observed` column | review chunk |
| `technology_labels` | `data/assays.csv` → `readouts` | review chunk |
| `target.construct_start/end`, `tags`, `expression_host` | `data/paper_constructs.csv` **(new)** → residue-slider bands | new feature |
| `plate.supplier`, `substrate.supplier` (catalogue numbers) | `data/kits.csv` | review chunk (`kit_proposals.csv`) |
| `instrumentation.readers`, `liquid_handling` | `build_output/equipment_candidates.csv` | existing report, new input |
| `interference_and_counterscreens.codes` | `data/assays.csv` → `interference` | review chunk |
| `paper.pmid / doi / title / year / journal` | `data/citations.csv`, and `assays.csv.learn_more_url` | `scripts/fetch_citations.py` |
| `unmapped_terms` | `scripts/term_table.csv` | review chunk — **and a relabel** |
| `not_reported` | a new gap report | tells us what the literature never states |
| `pharmacology`, `stimulation` | `data/assays.csv` → `effect`, and the entry's `conditions` text | review chunk |

### The field that is worth the most

`conditions`. Today `data/condition_stats.csv` is regex-mined from ChEMBL and PubChem assay
*descriptions*, and a gene-level row needs ten campaigns that state the parameter (`COND_MIN_N`). The
result: of the first 425 proteins surveyed, **295 have no gene-level condition rows at all** and show
class-level buffer and pH instead. A methods section states those values plainly and in full. This is
the single biggest, most mechanical gain.

### The field we could not get any other way

`quality.z_prime` and `quality.dmso_tolerance_*`. The AGM work concluded Z′ cannot be mined from
PubChem and that DMSO had to be a curated default. Papers do carry both — but thinly: in the 75-paper
sample only **4 of 73 stated a true Z′** and **26 stated a final DMSO**. So treat this as slow
accumulation, not a fix. One measured Z′ for one assay in one paper is a data point, not a standard.

### The genuinely new capability

`target.construct_start/end`. The app already reasons about constructs: `classStatus()` scores the
overlap between the user's residue range and annotated domains, and drops assays whose domain falls
outside it. Paper constructs are what people *actually express* — "WRN helicase core (500–942)",
"CDK2 (1–298)", "BLM helicase core (642–1296)". Showing those as bands on the residue slider, beside
the UniProt domains, would answer a question the app cannot answer today: *what construct do people
use for this assay?* 35 of 73 sample papers stated a range, so coverage would be partial but real.

## The pipeline, in five stages

```
build_output/papers/*.json              (one per paper, schema-validated)
        |
  1. tools/papers_aggregate.py          group by (entry_id, scope, parameter); n_papers, median,
        |                               spread, and the PMIDs behind each value
        v
  build_output/papers/paper_conditions.generated.csv
  build_output/papers/paper_constructs.generated.csv
  build_output/papers/paper_quality.generated.csv
  build_output/papers/term_candidates.csv        (from unmapped_terms)
  build_output/papers/kit_candidates.csv         (from plate/substrate suppliers)
        |
  2. tools/papers_to_review.py          -> review/NN_papers_<topic>/
        |                                  entry_proposals.csv, kit_proposals.csv,
        |                                  term_proposals.csv, REVIEW.md, decision columns
        v
  3. OWNER SIGNS OFF, then: tools/apply_review.py review/NN_papers_<topic> --apply
        |
  4. py scripts/build_database.py --aggregate-only --merge
     py tools/build_db.py
        |
  5. docs/assay_db.json -> the site
```

Two of those stages are new code (`papers_aggregate.py`, `papers_to_review.py`). Stages 3 to 5 are the
existing machinery unchanged, which is the point: this plugs into the review workflow that chunks 01–07
already use, rather than inventing a second way to change curated data.

### Aggregation rules, so one paper cannot move a default

- A value reaches a review proposal only with **at least 3 papers agreeing**, and the proposal carries
  `n_papers`, the median and the full spread. A single paper becomes a *note*, not a default.
- Entry-level columns are per **function class**, not per protein. `data/assays.csv` has 110 entries
  keyed by class; a TAAR1 cAMP measurement informs the GPCR cAMP entry as one observation among many.
  Aggregate per (entry_id, class) before proposing, never per protein.
- `confidence: low` and `extraction_flags` containing `method_by_reference` or `values_from_figure`
  are excluded from aggregation and listed separately for the owner to look at.
- `units_ambiguous` values are never aggregated. The TAAR1 IP1 assay is the example: the paper prints
  "5 ml" into a 384-well plate, which must mean µl. That belongs in front of a human, not in a median.

## What changes on the site

1. **"Reported by campaigns" gains a sibling.** `condFor()` in `docs/index.html` already resolves
   protein-level first, then class-level. It gains a third source so the details panel can show
   *Reported in methods sections (n = 4 papers)* beside the campaign-mined values, each linking to its
   PMID. Showing both, labelled, is better than replacing: campaign values are many and shallow, paper
   values are few and precise, and a screener wants to see where they disagree.
2. **Readout lists get literature evidence.** `readout_suggestions.csv` already flags 2,487
   class/readout pairs seen in campaigns but absent from every playbook entry. Paper
   `technology_labels` turn the strongest of those into review proposals with a citation attached.
3. **Quality target and DMSO cite measurements** where enough papers agree, instead of only the AGM
   default. The curated text stays; the measured range sits beside it.
4. **New: literature constructs on the residue slider.** Bands from `paper_constructs.csv`, visually
   distinct from the UniProt domain bands, snappable like them.

## Three things that will go wrong if ignored

### 1. Deduplication — and the map does not exist yet

A mined paper may already be a campaign. Campaign keys look like `doc:CHEMBL1155814` or `pc:686994`;
paper records carry a PMID and a DOI. **There is no DOI/PMID to ChEMBL-document map on disk** —
`build_cache/m1_data/chembl_docs/` does not exist.

The fix is already written but has never been run: `scripts/m2b_fetch_abstracts.py` fetches ChEMBL
document metadata and keeps `document_chembl_id`, `doi` and `pubmed_id`. Running it builds exactly the
join we need. **Treat that as a prerequisite** for anything that counts papers and campaigns together.

Until then the separate-tier design keeps us honest: double counting cannot happen if papers never
enter the campaign denominator.

Scale of the risk, for one protein: TAAR1's 46 campaigns are 22 `doc:CHEMBL*` and 24 `pc:*`, with
ChEMBL document IDs in the CHEMBL1155814–CHEMBL2203234 range, roughly 2010–2013. The ten papers we
mined are 2021–2026, so overlap there is probably nil. That is one protein and not evidence about the
other 5,139 — a heavily published target like EGFR, with 1,817 campaigns, will overlap much more.

### 2. Adding a vocabulary term forces a relabel

`unmapped_terms` is the growth path for `scripts/term_table.csv`, and it will fill up fast: the TAAR1
and WRN extractions alone produced 8 and 6 unmapped terms, including GloSensor, Transcreener
AMP2/GMP2 FP, the IP-One kit and ¹⁹F CPMG NMR. But a new technology term changes how every campaign is
labelled, so it needs the **full m2 relabel**, not just `--aggregate-only --merge`. That is the
expensive operation in this repo.

So: batch term additions into one review chunk and one relabel. Do not add them one paper at a time.

### 3. Paper values must not silently outrank curated ones

Everything in `data/assays.csv` was reviewed and signed off. A measured DMSO from three papers is
evidence to put *in front of the owner*, not a reason to overwrite a curated starting point. Hence
every route into a curated column goes through `entry_proposals.csv` with a `decision` column, exactly
like chunks 01–07, and `apply_review.py` already refuses to touch a cell the owner has changed since
the proposal was written.

## Suggested order of work

1. **`tools/papers_aggregate.py`** on the records we already have. Cheap, no network, and it shows
   what the yield per paper really is before more are mined.
2. **Run `scripts/m2b_fetch_abstracts.py`** to build the DOI/PMID join. Needed before any claim about
   papers and campaigns together, and it also unlocks the publication-year question left open from the
   readout work (whether radiometric leads for EGFR only because the literature spans 30 years).
3. **Mine 10–20 proteins from the top of `survey.csv`** with the agent, giving a seed set good enough
   to judge the local model against.
4. **`tools/papers_to_review.py`** and the first paper-derived review chunk — conditions only, because
   it is the highest volume and the least contentious.
5. **The site tier**, once there is enough in `paper_conditions.csv` to be worth rendering.
6. **Literature constructs on the slider** last: it is the most valuable but the only one that needs
   real UI work rather than a new row in a panel.
