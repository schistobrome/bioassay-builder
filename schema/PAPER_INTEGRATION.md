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
That reason stands on its own, and it is the whole reason: a paper is not a campaign, so it does not
belong in a count of campaigns. Double counting is a second, lesser worry, and as of 2026-11-18 it is
measured rather than assumed — 1 mined paper in 220 is a ChEMBL document behind a campaign for the same
protein (see "Deduplication" below). Small, but not zero, so the arithmetic objection is reduced rather
than removed. The semantic one is untouched.

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

### Built, 2026-11-18 — and what the plan above got wrong

`tools/papers_aggregate.py` and `tools/papers_to_review.py` now exist and were proven end to end on
batches 0–3 (408 records → `review/09_papers_batches_0_3`, then an accepted row through
`apply_review.py` and `build_db.py` against a copy of `data/`). Three corrections to the plan above:

1. **The counting tier was already built.** `scripts/literature_stats.py` turns the same records into
   `data/literature_stats.csv` — the "60% of 25 papers" number, per protein — and it works today
   (138 proposals from 220 extractions). So `papers_aggregate.py` does *not* duplicate it. It covers
   the other half: changes to `data/assays.csv` itself, which `literature_stats.py` cannot make.
2. **The input is `crawl/evidence.csv`, not the raw JSONs.** `scripts/paper_crawl.py log` has already
   resolved each record's protein classes, matched playbook entries and scope. Re-deriving that from
   `build_output/papers/*.json` would mean re-implementing the matching rule and risking it drifting.
3. **A new entry cannot be proposed without its ratings.** `build_db.py` refuses to build
   `docs/assay_db.json` *at all* if `screening_suitability` or `biological_relevance` is blank, so a
   single accepted row would take down the whole site build one step after the decision that caused
   it. `apply_review.py` now refuses such a row itself (`check_entry`), mirroring `build_db`'s
   vocabulary for `role`, `assay_category` and `scope` too. Papers cannot settle a rating — that is an
   expert judgement about whether an assay is worth running, not whether it is used — so the generator
   leaves both blank deliberately and the row is blocked until a person fills them in. This is the
   intended gate, not a defect.

The three candidate kinds the aggregator emits, in increasing cost to the owner: a **readout** is one
cell on an existing row; a **new entry** is a row that needs rating; a **term** needs the full m2
relabel and is not an `apply_review` operation at all.

### The `system` column, 2026-11-18

`data/assays.csv` gained a `system` column: the biological preparation an assay runs in — a Ba/F3
cytokine-independence line, HepaRG, an individual-donor microsome panel, an iPSC liver organoid. It is
distinct from `scope`, which says only cell-free or cell-based, and from `reagents`, which lists what
goes in the well. It is in `apply_review.py`'s `NEW_COLUMNS`, so the first accepted `system` proposal
creates it; `build_db.py` publishes it as `systems`; `docs/index.html` renders a **Systems used** line.

**Nothing can fill it yet, and the reason matters.** All 17 system and counter-screen values the crawl
found rest on **exactly one paper each**. That is not a shortage of papers but the shape of the
evidence: a methods section names its preparation once, and two papers on the same target rarely pick
the same one. Corroboration needs more papers *of the same class*, not more papers overall. This is the
strongest single argument for choosing later batches by empty class rather than by campaign count.

The literature *counting* tier needed no schema work at all: `data/usage_stats.csv` already has
`n_literature_total` and `n_literature_using`, `build_db.py` already folds them into the weighted
score at the `config.json` weight, and `docs/index.html` already renders "12% in papers" with an
outlined bar. Only the content was missing.

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

### 1. Deduplication — measured 2026-11-18: the overlap is real but rare

A mined paper may already be a campaign. Campaign keys look like `doc:CHEMBL1155814` or `pc:686994`;
paper records carry a PMID and a DOI, so the question needs a DOI/PMID to ChEMBL-document map. That map
now exists: `scripts/m2b_fetch_abstracts.py --crawled` has cached **all 12,867** documents behind the 80
crawled genes' campaigns (ChEMBL returns intermittent 500s, so it took three passes; each re-run
retries only the gaps and the third reported `1250 saved, 0 failed`). `tools/paper_campaign_overlap.py`
does the join.

**Result at complete coverage: 1 of 220 mined papers is a ChEMBL document behind a campaign for the
same protein.** That paper is PMC11215726 / PMID 38843875, `10.1021/acs.jmedchem.4c00197`, mined for
JAK2, curated by ChEMBL as CHEMBL5500392 with campaigns for JAK1, JAK2, JAK3 and TYK2.

The rate is low, not zero, and the distinction was earned the hard way: at 9,200 of the 12,867
documents the count was 0, and the next 2,417 produced the match (the final 1,250 added none). **Do not
read a zero from a partial cache as an absence of overlap.** The tool prints identifier coverage on
both sides, so a real absence can be told from a join with nothing to work with — but coverage of the
*document set* is the limit that actually bit, and that one is not visible from inside the join.

What this number does and does not cover: the cache holds the documents behind the **crawled genes'**
campaigns. Mined papers are mined for crawled genes, so the `same_gene` count is complete. The
`other_gene` count is not — a mined paper could sit behind a campaign for a gene outside the crawl set
and would show as `no_match`. Closing that needs `--all` (47,583 documents), and it matters only if the
tiers are ever combined, since an `other_gene` hit is not a double count for the protein we mined.

Why it is rare, with the mechanism stated correctly:

- **Years.** 134 of the 220 mined papers are 2025–2026. ChEMBL has 264 documents from 2025 and none
  from 2026 — curation lags publication by years. This half shrinks as ChEMBL catches up.
- **Journals.** ChEMBL's documents are medicinal-chemistry SAR papers: J Med Chem (3,152), Bioorg Med
  Chem Lett (2,618), Eur J Med Chem (1,029), Bioorg Med Chem, ACS Med Chem Lett. The crawl draws on the
  PMC open-access subset, which is biology and methods: Nature Communications, Science Advances, PLoS
  One, Scientific Reports, iScience. The corpora *mostly* miss each other — but the one match is itself
  a **J Med Chem** paper, because some ACS articles are open access and reach PMC. So the mechanism is
  "PMC open access rarely carries med-chem journals", not "never". Treat the overlap as a low rate to
  be measured, not a structural impossibility.

This still corrects the guess that stood here before, that a heavily published target like EGFR "will
overlap much more". At 1 in 220 it does not, and JAK2 — 1,000+ campaigns — is the one that hit.

Keep the tiers separate. Re-run `tools/paper_campaign_overlap.py` after any batch, and especially after
any change in where papers come from: a crawl that reaches closed-access med-chem papers would land
squarely in ChEMBL's corpus and the rate would climb.

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

**Status 2026-11-18: steps 1, 2 and 4 are done; step 3 is under way.** Steps 5 and 6 are open. The
crawl covered 80 targets across 4 batches, which is 15 of the 37 classes `data/assays.csv` is keyed
by — so step 3 is chosen **by empty class**, not from the top of `survey.csv`: another kinase informs
the same `K1` row as the first one did.

**Two of the 37 classes are not crawlable, and are not gaps.** No protein in `survey.csv` (5,852 rows)
carries either label, by design:

- `generic` — "Any protein (biophysical and interference checks)", noted "Always applies". It is
  applied to every protein, so no protein is labelled with it.
- `degrader_modality` — "Targeted protein degradation project", noted "Only included when you pass
  `--modality degrader`". It is selected by modality, not by what the protein is.

So the real denominator is **35 crawlable classes**, and the 22 "empty" ones are 20 gaps plus these
two. Counting them as gaps overstates the work and invites a search for proteins that cannot exist.

**Batch 4 (2026-11-18)** takes one target from each of the 20 reachable empty classes, which moves
class coverage from 15/35 to 35/35. Breadth first is deliberate: it measures yield per class before
committing the remaining slots. The plan beyond it is roughly three targets per class, because the
`system` evidence showed corroboration needs several papers **of the same class** — so batches 5 and 6
add depth to the classes batch 4 shows to be productive, rather than more breadth.

Targets were chosen for whether the class's characteristic assay will appear in the papers, not by hit
count, and the top-ranked candidate was often wrong for that: `SHH` ranks second in `nuclease` and is
not a nuclease, `COIL` ranks second in `chromatin_reader` and is not a reader, `HMGB1` and `FOXO1` top
`autophagy` without being autophagy machinery, and `NLRP3` tops `atpase` with 3,984 hits but its ATPase
activity is rarely the thing assayed. The picks are ABCB1, PARP1, HSPA5, SQSTM1, CBX8, BLM, USP7,
KDM1A, BRD4, EZH2, KAT8, CCND1, SPHK1, ESR1, APEX1, DNMT1, PDE5A, SNCA, KRAS, TOP1.

1. ~~**`tools/papers_aggregate.py`** on the records we already have.~~ Done. Cheap, no network, and it
   shows what the yield per paper really is before more are mined.
2. ~~**Run `scripts/m2b_fetch_abstracts.py`** to build the DOI/PMID join.~~ Done 2026-11-18. The
   dedup answer is above: 1 overlap in 220, at complete coverage. The join also unlocks the
   publication-year question left open from the readout work (whether radiometric leads for EGFR only
   because the literature spans 30 years) — the cache now carries a year per document.

   *2026-11-18:* the script's input path was stale — it looked for `scripts/m1_data`, while the cache
   has been at `build_cache/m1_data` since it stopped being a two-protein prototype, and it defaulted
   to the genes `WRN, EGFR`. It now takes `--crawled` (every target in `crawl/targets.csv`, ~12,900
   documents) or `--all` (every gene in the cache, ~47,600). `tools/paper_campaign_overlap.py` then
   answers the actual question and writes `build_output/papers/campaign_overlap.csv`: per mined paper,
   whether it is already a campaign for the same protein (`same_gene`, a real double count waiting to
   happen), for a different one (`other_gene`, not a double count — papers assay several targets), or
   absent.
3. **Mine 10–20 proteins from the top of `survey.csv`** with the agent, giving a seed set good enough
   to judge the local model against.
4. ~~**`tools/papers_to_review.py`** and the first paper-derived review chunk~~ Done, as
   `review/09_papers_batches_0_3`: readouts and new entries rather than conditions, because the crawl
   records carry function and technology labels but not `conditions.*`. Conditions remain the highest
   volume prize and still need either a crawl schema that captures them or the raw-JSON path.
5. **The site tier**, once there is enough in `paper_conditions.csv` to be worth rendering.
6. **Literature constructs on the slider** last: it is the most valuable but the only one that needs
   real UI work rather than a new row in a panel.
