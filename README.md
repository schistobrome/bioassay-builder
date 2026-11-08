# BioAssay Builder

App that scientists can use to simplify strategy for choosing which assays (biochemistry, biophysics, functional) to develop for a target protein and desired effect (e.g. inhibition, activation) that sources data from UniProt and PubChem Bioassay.

> [!NOTE]
> Work in progress! Project using Claude and public data to make assay selection easy and clear


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

## What "how standard" means
For one gene or function class: of the unique HTS campaigns whose assay type could be determined, the share that used each assay.
Shares do not add to 100% because one campaign can use several assays. The score mixes campaign share and literature share
using the weights in `data/config.json`. Rows with fewer than 10 classifiable campaigns are flagged low confidence.
Public assay text is often too sparse to determine technology, so always show "based on k of N campaigns".