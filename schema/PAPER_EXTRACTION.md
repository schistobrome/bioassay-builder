# Mining assay detail out of open-access papers

Two deliverables, both built against the 75 sample PDFs in `data/public/pdf` and tested against the live
Europe PMC API on 2026-11-15:

| File | What it is |
|---|---|
| `schema/paper_extraction.schema.json` | The extraction format. One record per paper, with an `assays` array. |
| `schema/paper_extraction.example.json` | A real filled record for an **inhibitor** paper (PMC12770528, three assays). |
| `schema/paper_extraction.example-activator.json` | A real filled record for an **activator / PAM** paper (the GLP-1R duplexed screen, two assays). |
| `data/paper_queries.json` | The search configuration: query shape, per-class assay terms, exclusions, re-ranking recipe. |
| `scripts/fetch_papers.py` | The runnable pipeline: search, fetch, score, split into assay blocks, and optionally extract with a local LLM. |

Both files validate: the schema is valid JSON Schema 2020-12, the example validates against it, and its
vocabularies were checked against `scripts/term_table.csv`, `data/assays.csv` and `data/function_classes.csv`.

## What the 75 papers actually contain

I extracted text from 73 of the 75 (two failed with a pypdf stream-token limit) and counted what each one states.

**Nearly always there:**

| Detail | Papers (of 73) |
|---|---|
| Detergent | 72 |
| A named detection technology | 66 |
| Analysis software | 62 |
| Incubation time / temperature | 61 |
| A potency output (IC50 / EC50 / Kd / Ki) | 61 |
| Buffering agent | 59 |
| Protein vendor | 53 |
| Named instrument | 52 |
| Plate format | 50 |

**Often missing, and these are the ones the app most wants:**

| Detail | Papers (of 73) |
|---|---|
| Protein construct residue range | 35 |
| Explicit final DMSO | 26 |
| Z′ or signal-to-background | 26 (only 4 state a true Z′) |
| Reagent catalogue number | 19 |
| Assay volume | 15 |
| Plate vendor + catalogue number | 11 |
| Enzyme concentration | 9 |

So the realistic expectation is: readout, buffer, plate format and potency come out of most papers; construct
boundaries and DMSO out of roughly a third; Z′ out of almost none. That matches what the AGM work already
concluded — a curated default is the only workable route for DMSO, and Z′ cannot be mined at scale.

**One paper describes many assays.** Method-heading counts ranged from 1 to 10. The richest
(`s41586-024-07318-y`, WRN) had 10: helicase unwinding, ATPase, HTRF, protein expression, GSH reactivity,
high-content nuclear count and more. The schema therefore has an `assays` array, not one assay per paper.

## Why the extraction needs an LLM and not regex

I tried the obvious regex for construct ranges. On these papers `(\d+[-–]\d+)` matched `(500-942)` correctly for
the WRN helicase core, but it also matched `(2015–2018)`, `(12–17)`, `(42–44)` and dozens of other citation
years and page ranges. Pulling a residue range needs the cue word nearby and the protein it belongs to, which
is a judgment, not a pattern. The same goes for:

- **Stock vs final concentration.** "compounds dissolved in 100% DMSO" is the stock; `dmso_pct` wants the 1% in the well.
- **Which readout belongs to which assay.** A paper with an HTRF kinase assay and a KINOMEscan panel must not
  credit KINOMEscan to the kinase assay. This is the same trap that distorted the readout counts in the app.
- **Artefact corrections that look like data.** The SPR section of PMC12770528 builds a DMSO curve from 0% to 2%.
  That is a refractive-index correction, not a DMSO tolerance test. The example record leaves
  `quality.dmso_tolerance_tested_pct` empty on purpose and says why.

## Activators and inhibitors

My first draft of this was quietly inhibitor-biased, and it needed fixing. Counting the activation
vocabulary across the sample: **11 of the 73 papers are substantially activation studies** — the GLP-1R PAM
and agonist screens, GPCR biased agonism, CDK2 and STAT6 activators. `potentiator` appeared in 9 papers and
`positive allosteric modulator` in 6, and **none of those would have matched any inhibition search term.**

What the schema now carries for activation, beyond `effect_direction: activation | both`:

| Field | Why an inhibition-only schema loses it |
|---|---|
| `pharmacology` | "activation" covers designs that are not interchangeable. A direct agonist screen runs on unstimulated cells; a PAM screen only works against a sub-maximal agonist challenge. |
| `stimulation` | The reference agonist and its level (EC20 / EC50 / EC80 / ECmax / unstimulated). **An EC50 from a potentiator screen is meaningless without it** — potency at EC20 agonist is not potency at EC80. |
| `protocol_additions` | A triple-add GPCR screen is one assay with three reads, not three assays. |
| `controls.basal`, `controls.max_stimulation` | An activation assay's floor is the unstimulated well and its ceiling is ECmax agonist or forskolin. Those are not the positive/negative controls of an inhibition assay. |
| `quality.basal_signal`, `quality.fold_response` | High basal leaves no room to go up: the activation equivalent of a poor signal window. |
| `outputs` | Added `EC20`, `EC80`, `Emax`, `percent_activation`, `fold_induction`, `pEC50`, `hill_slope`, `Kb`, `DC50`, `Dmax`. |
| `screening_campaign.hit_categories` | One screen can sort compounds into agonist / potentiator / antagonist. |

The search terms gained `agonist`, `partial agonist`, `inverse agonist`, `positive allosteric modulator`,
`allosteric modulator`, `potentiator`, `enzyme activator`, `EC50`, `Emax`, `fold induction` and
`maximal response`; the GPCR block gained `EC20`, `EC80`, `triple-add`, `three-phase assay` and
`biased agonism`; and the scorer gained an `activation_design` signal. There is deliberately **no exclusion**
for agonist or activator work.

`schema/paper_extraction.example-activator.json` is the proof: it is the real GLP-1R duplexed PAM screen, with
the EC20 challenge, the three-add protocol, the vehicle / EC20 / EC80 / ECmax control row, Z′ = 0.68, 175,000
compounds and 98 hits sorted into agonist, potentiator and antagonist.

## The schema, and the three decisions behind it

**1. It speaks the app's existing vocabulary.** Nothing new is invented:

| Schema field | Comes from |
|---|---|
| `function_labels` | the 55 function terms in `scripts/term_table.csv` |
| `technology_labels` | the 32 technology terms in the same file |
| `assay_category`, `role`, `scope` | `data/assays.csv` |
| `conditions.*` | `COND_NUMERIC` / `COND_CATEGORICAL` in `scripts/build_database.py` — `ph`, `incubation_min`, `dmso_pct`, `mgcl2_mm`, `nacl_mm`, `atp_um`, `bsa_pct`, `volume_ul`, `temperature`, `buffer`, `detergent`, `reducing_agent` |
| `quality.qc_metric_code`, `interference_and_counterscreens.codes` | the AGM columns added in review chunk 07 |
| `instrumentation.instrument_classes` | the `INSTL` keys in `docs/index.html` |

Because the condition keys match exactly, extracted values merge into `condition_stats.csv` with no translation
layer. Where a paper's units vary too much for a fixed-unit key — enzyme at 250 pM, substrate at 50 µM,
compound stock at 100 mM — the schema uses a `{value, unit, qualifier, verbatim}` object and keeps the paper's
own unit rather than converting.

**2. Absence is recorded, not implied.** Every assay has a `not_reported` array of dotted field paths. A null with
no entry there means the model missed it; a null with an entry means the paper does not say. Without this you
cannot tell a gap in the literature from a gap in the extraction, and the DMSO and Z′ counts above are exactly
the kind of question that needs the difference.

**3. Every assay carries verbatim quotes.** `evidence` is required and holds at least one exact quote tied to a
field path. This is what makes the output reviewable the same way the `review/` chunks are: a scientist can check
any number without reopening the PDF. `extraction_flags` carries the caveats, and the most important is
`method_by_reference` — when a paper says "as described previously", the conditions belong to the cited paper,
not this one, and silently attributing them would corrupt the aggregate.

Terms the model cannot map go in `unmapped_terms` rather than being forced into the nearest enum. In the worked
example, "19F CPMG relaxation NMR" matches none of the 32 technology labels, so it lands there with its quote.
That list is the curation queue — it is how the vocabulary grows, in the same sign-off pattern as the review chunks.

### Running it with a small model

Do not hand a small model the whole paper and the whole schema. Three passes:

1. **Segment** — no model needed. Take `<sec><title>` from the Europe PMC full text XML. On the worked example
   that gave 16 methods subsections, the ATPase block being 1,922 characters.
2. **Extract** — one subsection at a time, asking only for one assay object. Give the model the enum lists for
   `function_labels` and `technology_labels` in the prompt, and tell it to leave anything unstated null.
3. **Map and merge** — reconcile `target` across assays (the purification section usually holds the construct and
   tags for assays described later), assign `assay_id`s, and collect `unmapped_terms`.

Two things to prompt for explicitly, because they are where small models fail: final rather than stock
concentrations, and never copying a value from a neighbouring assay's buffer.

## The search queries

### Europe PMC, not PubChem

PubChem has no full text, and PubMed serves abstracts. Europe PMC is the only one of the three with an API that
returns open-access full text, and it supports section-scoped fields, which is what makes these queries selective.
Both endpoints are in `data/paper_queries.json`.

### The query shape

```
({target names in TITLE_ABS}) AND ({assay terms in METHODS}) AND (OPEN_ACCESS:y AND IN_EPMC:y) {exclusions}
```

The two scoping choices carry all the selectivity, measured on WRN:

| Query | Hits | What was on top |
|---|---|---|
| `WRN AND inhibitor`, free text | 6,933 | DNA-damage reviews from 2009–2014 |
| target free-text + assay terms | 869 | still reviews |
| **target in TITLE_ABS + assay terms in METHODS** | **22** | the four WRN inhibitor-discovery papers |

Putting the target in `TITLE_ABS` makes the paper be *about* the protein; putting the assay terms in `METHODS`
makes the match a real protocol rather than a passing mention. Verified live:

- **WRN** → 18–22 hits, led by the on-target resistance paper, the structural/conformational paper,
  the VVD-214 covalent allosteric inhibitor paper and the fragment-based lead discovery paper.
- **EGFR** → 1,239 hits including "A Fluorescence-Based High-Throughput Screening Method for EGFR Inhibitors".
- **BRD4** → 182 hits including "Fragment-Based Discovery of BRD4-BD2 Selective Binders: Biophysical Validation".
- **LRRK2** → 98 hits including docking and HTS papers.

**Do not sort by citation count.** `sort=CITED desc` put 2009–2014 reviews above every assay paper, because old
reviews accumulate citations and mention every term. Leave sort empty for relevance.

### Per-class assay terms

`assay_terms_by_class` covers all 37 `class_id` values in `data/function_classes.csv`, so a script can look up the
classes the app already assigned to a protein and union the matching blocks with `assay_terms_core`. A kinase gets
ADP-Glo, Kinase-Glo, mobility shift and ATP Km; a GPCR gets cAMP, calcium mobilisation, β-arrestin, GloSensor and
GTPγS; an ion channel gets patch clamp, thallium flux and rubidium efflux. All five blocks I spot-checked returned
on-target assay papers.

### Relevance alone will not give you a top 10

Even with the scoped query, relevance ordering mixes in target-biology papers — for CDK2 it put "Benzoylaconine
suppresses colorectal cancer progression" above the allosteric-inhibitor assay paper. So the config specifies:
fetch 50 candidates, score each full text for assay-detail density, keep the best 10, and only spend the model on
papers scoring 4 or more.

The scorer gives a point each for a stated plate format, detection technology, buffer composition, potency output,
protein construct and named instrument, two points each for a quality metric or an explicit final DMSO, and a
penalty for biologics-only papers. Checked against the 75-paper sample: **22 papers scored 5 or 6, and all 22 were
small-molecule papers.** Mean score was 3.6 for small-molecule papers against 2.0 for biologics.

### Expect nothing for antibody targets

Of the 75 samples, 20 were antibody, CAR-T or radioligand papers — CD19, CD20, BCMA, PSMA, CD3 — and **not one
reached 5 of 6 assay-detail signals.** The `exclusions` block drops that literature at query time, and a CD20
query returns 43 mostly irrelevant hits rather than noise. An empty result for a biologics target is the correct
answer, not a bug. Targets where the sample was richest were WRN, EGFR, BRD4, CDK2, LRRK2, GLP-1R and STAT6.

### Licences

Europe PMC returns a `license` field: the WRN papers tested were `cc by`, `cc by-nc` and `cc by-nc-nd`. Extracting
facts is fine for all of them; redistributing the text of an ND paper is not. The schema keeps `paper.license` so
that stays visible, and the extraction stores quotes rather than whole articles.

## Running it

`scripts/fetch_papers.py` does all of it. Standard library only for the first four stages.

```powershell
# see the queries, no network
py scripts\fetch_papers.py --genes WRN,EGFR,GLP1R --dry-run

# search, fetch, score, split into assay blocks
py scripts\fetch_papers.py --genes WRN,EGFR,GLP1R --keep 10

# every protein the app knows, resumable
py scripts\fetch_papers.py --all --keep 10 --max-minutes 120
```

Then extract with a local model. Ollama is the default backend:

```powershell
winget install --id Ollama.Ollama --source winget --accept-package-agreements --accept-source-agreements
# open a NEW PowerShell window so PATH is picked up, then:
ollama pull qwen2.5:7b-instruct
pip install jsonschema

py scripts\fetch_papers.py --genes WRN --extract --model qwen2.5:7b-instruct
```

On Windows the installer starts the server itself (tray app / service), so **do not run `ollama serve`** —
it fails with "address already in use". Check it is up with `ollama list`, or
`Invoke-WebRequest http://localhost:11434/api/tags`.

**Pick the model to fit your VRAM**, since the whole point is that the server holds it resident:

| VRAM | Model | Roughly |
|---|---|---|
| 8 GB | `qwen2.5:7b-instruct` | ~4.7 GB at Q4_K_M, leaves room for a 16k context |
| 12 GB | `qwen2.5:14b-instruct` | ~9 GB at Q4_K_M |
| 24 GB+ | `qwen2.5:32b-instruct` | best at the stock-vs-final and don't-borrow-the-neighbouring-buffer rules |

A 14B model on an 8 GB card partly offloads to CPU and runs several times slower; it still works, so it is
worth one timed paper before committing to a long run. If a block overruns the context, lower
`--block-chars` (default 12000) rather than raising `--num-ctx` past what the card holds.

Any OpenAI-compatible server (llama.cpp `--server`, vLLM, LM Studio) works too:

```powershell
py scripts\fetch_papers.py --genes WRN --extract --backend openai `
   --host http://localhost:8000 --model <model-name>
```

Each block goes to the model on its own with `temperature 0` and JSON mode, and the finished record is
validated against the schema before it is written. Output is one file per paper in `build_output/papers/`,
plus `index.json`; caches are in `build_cache/papers/{search,fulltext,blocks}`, so reruns cost nothing.

### Or let the Claude Code agent fill it

`--backend claude` is the third route, for the handful of papers where extraction quality matters more
than throughput. The script cannot call Claude Code — Claude Code is what runs the script — so this
inverts the direction: the script does search, fetch, score, split, merge and validate, and the agent
does the reading and judging.

```powershell
# 1. prepare the queue
py scripts\fetch_papers.py --genes WRN --keep 2 --extract --backend claude
```

That writes `build_output/papers/agent_queue/`:

```
agent_queue/TASKS.md              instructions + a table of every block to do
agent_queue/queue.json            the same, machine-readable
agent_queue/<pmcid>/paper.json    pmid / doi / title / journal / license / score, so the agent
                                  never retypes paper metadata
agent_queue/<pmcid>/A1.json       <- the agent writes these, one assay object each
```

```powershell
# 2. in Claude Code: "Work the queue in build_output/papers/agent_queue/TASKS.md"
# 3. assemble, hoist unmapped_terms, validate
py scripts\fetch_papers.py --agent-merge

# 4. check any time
py scripts\fetch_papers.py --validate build_output/papers
```

`TASKS.md` carries the same eight rules the local-model prompt carries, so the two routes produce
comparable records. A block whose JSON already exists is marked done and skipped, so it is resumable:
do a few, merge, come back. `--agent-merge` is idempotent, and refuses a file that looks like a whole
paper record instead of one assay.

**Why one file per assay** rather than one per paper: a block that goes wrong loses one assay, not the
paper; the agent can stop and resume anywhere; and `assay_id` is derived from the filename, so the
ordering cannot drift from the block numbering.

**Worked for real on PMC10358344** (the 2-sulfonyl/sulfonamide WRN covalent inhibitor paper): three
blocks, agent-filled, merged clean against the schema.

| assay | technologies | flags |
|---|---|---|
| A1 Protein Expression and Purification | — (not an assay) | `values_from_supplementary` |
| A2 Multiplexed WRN Activity Assay | Fluorogenic / FRET substrate + Fluorescence polarization | `multiple_assays_merged_in_source`, `units_ambiguous` |
| A3 TAMRA-Dequenching DNA Unwinding | Fluorogenic / FRET substrate | `method_by_reference` |

Three things in that record are worth looking at, because they are what the small-model route most often
gets wrong:

- **A3 is flagged `method_by_reference`** and marked `confidence: medium`, because its buffer, plate and
  reader settings exist only in A2 ("the same as described above"). The values are carried over and the
  flag says so, rather than being silently attributed to A3 as if measured there.
- **`0.1 mg/mL BSA` did not become `bsa_pct`.** That fixed-unit key wants a percentage, and mg/mL is not
  convertible without a density assumption, so `bsa_pct` is null, the value sits in
  `conditions.other_components`, and the mismatch is recorded in `unmapped_terms`.
- **A1 records no assay at all** — it is reagent preparation. It still carries the ZnSO4 supplementation
  and the expression hosts, which is real construct-development detail, with `function_labels` and
  `technology_labels` left empty rather than guessed.

**When to use which route**

| | local model | Claude Code agent |
|---|---|---|
| Cost per paper | near zero | session tokens; the expensive path |
| Throughput | thousands of blocks | tens |
| Judgment calls (stock vs final, borrowed buffers, method-by-reference) | needs grammar constraint and tuning | handled |
| Good for | the `--all` sweep | a seed set, and a gold standard to tune the local prompt against |

The sensible pattern is both: have the agent do 10–20 papers properly, then use those records as the
reference for judging whether the local model is good enough to turn loose on the sweep.

**Model size.** The per-block prompt is roughly 2–4k tokens of instructions and vocabulary plus up to 12k
characters of methods text, so give the model at least a 16k context (`--num-ctx`, the default). A 7B instruct
model will follow the shape; the fields that reward a 14B or larger are the stock-versus-final distinction and
not borrowing a neighbouring assay's buffer.

### What I tested, and what I did not

Tested live on 2026-11-15:

- **Search, fetch, score, split**, for real. WRN returned 94 hits; the top six by score were the covalent
  inhibitor paper (10), the on-target resistance paper, the structural paper, the chromatin-degradation paper,
  the fragment-based discovery paper and the Picco WRN inhibitor paper (8 each). GLP1R returned 948 hits and
  surfaced "Novel Small Molecule GLP-1R Agonists" — the activator path working. 12 papers gave 71 assay blocks.
- **Block splitting.** PMC10358344 split into 9 blocks: protein expression, multiplexed WRN activity,
  TAMRA-dequenching unwinding, ADP-Glo ATPase, PK/LDH coupled, two exonuclease assays, ITC. The ADP-Glo block
  was 1,276 characters with the full buffer, plate catalogue number and ATP concentration intact.
- **The extract wiring**, against a mock server: prompt assembly, the call, JSON parsing, `assay_id` assignment,
  hoisting `unmapped_terms` up to the record, and schema validation. The finished record validated.

Not tested: **a real local LLM.** No Ollama or OpenAI-compatible server is installed on this machine, so the
extraction quality is unmeasured. Treat `--extract` as wired-and-validated, not as tuned. Start with one paper
and read the JSON yourself.

Two limitations worth knowing:

- **The scorer is less sharp on GPCRs than on enzymes.** For GLP1R it gave 8 to a humanised-mouse paper and 6 to
  an in-vivo profiling paper, because GPCR pharmacology papers carry assay vocabulary even when the assay detail
  is thin. WRN, where the scores were 8–10 for genuine assay papers, is the better case. Expect to raise
  `--min-score` for receptor targets.
- **The `--all` run is large.** 5,852 proteins at a 50-paper pool is a lot of API calls, so use `--max-minutes`
  and let it resume. Start with 20 proteins.
