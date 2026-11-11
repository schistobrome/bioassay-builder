# BioAssay Builder

App that scientists can use to simplify strategy for choosing which assays (biochemistry, biophysics, functional) to develop for a target protein and desired effect (e.g. inhibition, activation) that sources data from UniProt and PubChem Bioassay.

> [!NOTE]
> Work in progress! Project using Claude and public data to make assay selection easy and clear

[https://schistobrome.github.io/bioassay-builder/](https://schistobrome.github.io/bioassay-builder/)

The data lives in plain tables (CSV files) in `data/`. A script checks them and builds one file, `docs/assay_db.json`, which the website reads.

| File | Item |
|---|---|
| `data/assays.csv` | one assay type (a generic method) |
| `data/kits.csv` | one commercial kit or service |
| `data/usage_stats.csv` | how often assay is used in one scope (gene or function class) |
| `data/citations.csv` | a paper or study supporting assay or mechanism |
| `data/function_classes.csv` | a protein functional domain or class |
| `data/config.json` | settings (weights, thresholds) |

Lists inside a cell are separated by `;`. Levels are exactly: very_low, low, moderate, high, very_high.

## Helpers (`scripts`)
- `scripts/find_literature.py "BRD4 AlphaScreen"` finds papers (Europe PMC) and writes candidates to consider adding to data files
- `scripts/fetch_citations.py` fills titles, years and journals for the DOIs you put in `assays.csv`
- `scripts/usage_from_m2.py EGFR --m2-dir <folder>` turns labeled campaigns sourced from `m2_label.py` into usage-share rows
- `scripts/m*` staged sourcing and munging of unstructured/structured data from UniProt, CHEMBL, and PubChem BioAssay

## Running the data pipeline (Windows PowerShell)

Run the steps in this order from the repo folder (`cd "C:\Users\xanox\Documents\GitHub\bioassay-builder"`).
Every step can be stopped (Ctrl+C) and rerun: finished work is skipped. `build_cache/`, `build_output/` and
`data/public/` are large and can be regenerated, so never commit them.

| Step | Command | Reads | Writes |
|---|---|---|---|
| 1. PubChem bulk (optional, faster than the API) | `py scripts\download_bulk.py --url https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/ --dest data\public\pubchem\Extras --only bioassays.tsv.gz,RevokedAids.gz` then `py scripts\pubchem_plan.py`, the download command it prints, then `py scripts\pubchem_local.py` | `data\public\pubchem\` | `build_cache\pubchem_bulk\` |
| 2. Fetch, label, summarize | `py scripts\build_database.py --max-minutes 240` (add `--pubchem-local` after step 1) | UniProt, ChEMBL, PubChem (internet or step 1), `data\function_classes.csv`, `scripts\term_table.csv` | `build_cache\state\`, `m1_data\`, `m2_data\`, `summary\` |
| 3. BindingDB (optional) | `py scripts\bindingdb_label.py --folder data\public\bindingdb` | BindingDB TSVs, `build_cache\m2_data\` | `build_cache\summary_bdb\` |
| 4. Plate formats and conditions | `py scripts\plate_formats.py` and `py scripts\conditions.py` | `build_cache\m1_data\` | `build_cache\plates\`, `build_cache\conditions\` |
| 5. Aggregate and merge | `py scripts\build_database.py --aggregate-only --merge --replace-m2` | everything in `build_cache\`, `data\assays.csv` | `build_output\*.csv`, `data\usage_stats.csv`, `data\condition_stats.csv` (a backup goes to `data\backup\`) |
| 6. Website data | `py tools\build_db.py` | `data\*.csv`, `data\config.json` | `docs\assay_db.json` |
| 7. Checks (optional) | `py scripts\equipment_scan.py`, `py tools\audit_gaps.py`, `py -m unittest discover -s tests` | `build_cache\`, `build_output\`, `data\` | `build_output\equipment_*.csv`, `build_output\audit\` |

`scripts\build_pipeline.py` is an older name for step 2 and translates old options.

**If something fails**
- *Network errors or "STOPPED: the connection kept failing":* rerun the same command later. Downloads are cached in small pieces, so nothing is lost.
- *Proteins "parked":* they were too big for the first pass. Rerun with `--heavy-minutes 90`. `py scripts\build_database.py --status` lists them, and so does `build_output\parked_proteins.csv`.
- *"quarantined" (failed 3 runs):* `py scripts\build_database.py --retry-quarantined`.
- *You edited `scripts\term_table.csv`:* the next run relabels without downloading.
- *`data\usage_stats.csv` looks wrong after a merge:* copy the newest file from `data\backup\` back over it.
- *`tools\build_db.py` says "Not built":* it lists each ERROR row. Fix that row in `data\` and rerun.

## What "how standard" means
For one gene or function class: of the unique HTS campaigns whose assay type could be determined, the share that used each assay.
Shares do not add to 100% because one campaign can use several assays. The score mixes campaign share and literature share
using the weights in `data/config.json`. Rows with fewer than 10 classifiable campaigns are flagged low confidence.
Public assay text is often too sparse to determine technology, so always show "based on k of N campaigns".