# Code and data audit: Assay Builder

Branch `audit-and-gaps`, 2026-11-10/11. All changes are uncommitted. Knowledge gaps are in [GAPS.md](GAPS.md).

## Verdict

The pipeline is careful and resumable, and the front end works end to end. A headless browser run of EGFR, GLP1R, KCNH2, ESR1 and BRD4 in both directions produced no script errors. The main problem is in the numbers, not the code. **The same patent study is often counted two or three times**: once as a PubChem deposit, once as a BindingDB row and once as a ChEMBL document. That inflates campaign counts by 2–14% per class and shifts 122 of the usage shares on the site by more than 5 points. This branch fixes it in code. To get corrected numbers on the site, you must run three local commands (see "How to run"); nothing downloads. The second problem is presentation. The headline "how standard" share measures how often a *function* is studied, so sibling assays (patch clamp vs thallium flux, FRET vs gel) always show the same number, and 71% of proteins quietly use class-wide figures. Both need your decision, because the clean fix touches `assays.csv`. I also fixed several smaller bugs: the page had no mobile viewport, stale rows were never pruned, the equipment scan undercounted, and the readout suggestions raised false alarms. Twenty-one offline tests cover the code I touched, and all pass. The site JSON is 12.6 MB (0.5 MB as served compressed), under the 15 MB budget.

## Issues found

Likelihood is how often the problem actually occurs in the current data; impact is what it does to what a scientist sees.

| # | Severity | File | Issue (likelihood / impact) | Status |
|---|---|---|---|---|
| 1 | **High** | `scripts/build_database.py`, `scripts/bindingdb_label.py` | **The same study is counted more than once.** PubChem patent keys keep the kind code (`patent:US9023767B2`) and BindingDB keys do not (`patent:US9023767`), giving 3,070 duplicate patents. ChEMBL source 37 republishes BindingDB patents under a ChEMBL document id, giving 7,104 one-to-one duplicates by identical assay text (9,899 before the one-to-one filter). *Likely: about half of all BindingDB campaigns. Impact: class totals inflated 2–14%; 122 shares off by more than 5 points, 28 by more than 10.* | **Fixed.** Shared `canonical_campaign()` normalizes keys; BindingDB summaries now carry `aliases` (ChEMBL doc → patent). Plates and conditions are attached before re-keying. Tests added. **Data not yet updated: run the How-to-run steps.** |
| 2 | Medium | `build_database.py` `merge_into_data` | Automatic rows for a protein that now produces no rows were never removed. 83 `bd-20261109` rows (LGALS8, AGPS, SAP18 …) are still shown through the site's "all directions" fallback. *Certain; small.* | **Fixed**: summarized proteins count as regenerated, so their stale rows are pruned at the next `--merge`. Hand-made rows are never touched (test). |
| 3 | Medium | `build_database.py` readout suggestions | Generic entries (DSF, SPR, ITC) were not counted as covering any class, so `readout_suggestions.csv` flagged them as missing everywhere (about 3,400 campaign-pairs). *Certain; noise in a review file.* | **Fixed** + test. |
| 4 | Medium | `scripts/equipment_scan.py` | Only the first record of each campaign was scanned (often the short ChEMBL one-liner, not the PubChem protocol), and a campaign shared by two proteins counted only for the first protein. *Common; undercounts equipment.* | **Fixed** + test. |
| 5 | Low-medium | `scripts/equipment_scan.py` | Seed regexes matched ordinary words: `gemini`→geminin, `\bm5\b`→M5 muscarinic receptor, `mosquito`→mosquito cell lines, `harvest`→"cells were harvested", `omega`, `synergy`, `victor`, `prominence`, plain `hamilton`/`multichannel`, and `molecular devices` (which also makes FLIPR). *Common in the diagnostic file only; the site does not use it.* | **Fixed** + test. Rerun `equipment_scan.py`. |
| 6 | **High (phones)** | `docs/index.html` | No `<!doctype>`, `<meta charset>` or viewport tag, so phones rendered the desktop layout zoomed out and the browser ran in quirks mode. | **Fixed.** Verified headless at 390 px wide: no sideways scroll. |
| 7 | Low | `docs/index.html` | Header showed the literal placeholder "last updated YYYY-MM-DD" (no `reviewed_on` dates exist). | **Fixed**: says "no review date recorded yet · data built …". |
| 8 | Low | `docs/index.html` | Gene name put into a regex unescaped in `kitsFor` (names with `.`). Possible crash on "How standard" details when a hand-made row has blank counts (`level.replace` on null). | **Fixed.** |
| 9 | Low (a11y) | `docs/index.html` | Plan drawer had no Escape key or focus handling; lookup status was not announced to screen readers; dead variables in `condFor`. | **Fixed** (focus moves into and out of the drawer, Escape closes it, `aria-live` added). |
| 10 | Low (latent) | `build_database.py` `merge_effect` | A blank effect made a campaign "mixed", which counts it in both the inhibition and activation views. None seen in current summaries. | **Fixed** + test. |
| 11 | Medium (legacy) | `scripts/usage_from_m2.py` | Counted every assay row as a campaign (a paper with 12 assays counted 12 times). | **Fixed** + test. |
| 12 | Maintenance | `build_database.py`, `plate_formats.py`, `conditions.py`, `equipment_scan.py` | The "ChEMBL copy joins its PubChem original" logic was copied four times; the class-level merge duplicated `merge_campaign_lists`. | **Fixed**: one helper (`pubchem_originals` / `record_campaign`). Verified identical output on 49 real proteins for plates and conditions. |
| 13 | Low | `tools/build_db.py` | Entry-id set rebuilt for every condition row; no warning as the JSON grows. | **Fixed**: warns above 15 MB. Output verified byte-equivalent to the committed `assay_db.json` (except the date). |
| 14 | **High** | `build_database.py` usage rows | **Headline share is per function, not per format.** 21 of 72 entries always show the same headline number as a sibling (IC1/IC2/IC3, P1/P2, T1/T2, H2/H3 …). On hERG, patch clamp, thallium flux and membrane-potential dye all show "94%". *Certain; misleading.* | **Open: decision D2.** |
| 15 | Medium | `scripts/m2_label.py` | Effect is guessed from the whole text with "inhibit" checked first; it misses "potentiate", "positive allosteric", "induce" and "opener". About 2,200 GPCR, NR and ion-channel records whose *title* says agonist or potentiator are "unspecified", and unspecified counts toward inhibition. BindingDB already uses a better title-first rule (`bindingdb_local.effect_of`). *About 1–4% of records in those classes.* | **Open: D3** (changing it means relabeling). |
| 16 | Medium | `build_database.py` aggregate | A protein in two classes lends all its campaigns to both (an autophagy kinase's kinase campaigns appear in "autophagy" totals). Class shares are diluted, and class `n_campaigns` and the readout suggestions are noisy. | **Open: D4.** |
| 17 | Medium | `docs/index.html` `usageFor` | Below 5 classifiable campaigns the card silently switches to class-wide data (a small pill only); if no tier reaches 5 it shows the first tier anyway (for example "100% of 1"). 71% of proteins are below 5. | **Open: D5.** |
| 18 | Medium | `data/function_classes.csv` | The `kinase` keyword puts metabolic kinases and chemokines (PGK2, MVK, CCL3) into protein_kinase. GO `channel activity` puts aquaporins, claudins and UCP1 into ion_channel. | **Open: D6** (curated file, not edited). |
| 19 | Medium | `bindingdb_label.py` | BindingDB's own literature rows (`pmid:` keys, 4,227) may duplicate ChEMBL papers (`doc:` keys). This can't be measured without a ChEMBL document→PMID table; no word-for-word text matches were found. | **Open: D7.** |
| 20 | Low | repo | `build_output/` (including a 30 MB CSV) and `data/backup/` are tracked in git despite "never commit". There is a stray second ignore file named `gitignore` (no dot). | **Open: D1.** |
| 21 | Low | `docs/assay_db.json` | 12.6 MB raw / 0.5 MB gzip; 97% is usage rows that repeat the same `source` and `notes` text. Fine now; about 20% more rows would cross 15 MB. | **Open: D8.** |
| 22 | Low | `build_database.py` `conditions_in` | First stated value wins, so a detection or wash buffer listed first is taken as the assay buffer. `\d{2} c` can read "25 c…" fragments. Plausibility limits are applied. The site already says these are descriptive. | Pinned by tests; no change. |
| 23 | Low | `build_database.py` | Technology-row `stat_id` uses the first 12 letters of the label, so two labels could collide. No collisions today, and the site does not use `stat_id`. | Noted. |
| 24 | Low (a11y) | `docs/index.html` | Dimmed cards use 55% opacity, which can fall below text-contrast guidelines. | Noted (design choice). |
| 25 | Low | `scripts/pubchem_local.py` | Writes the same description into both `summary_text` and `full_description`. Harmless (labels use sets), but it doubles the text. | Noted. |

Checked and found correct: plate-format extraction ignores compound, source and stock plates (tests); `effect_subset` matches what the Definitions panel says; the merge never overwrites hand-made rows (test); the cell-based vs cell-free match between entries and campaigns works; every technology label in the data is classified by the instrument keyword map (only 3 playbook readout phrases are not: IC2 "rubidium efflux with atomic absorption", IC5 "binding", E9 "Competition binding panel").

**Output formats:** no column of `usage_stats.csv`, `condition_stats.csv` or `assay_db.json` changed. One addition: `build_cache/summary_bdb/<acc>.json` gains an optional `"aliases"` key. Older files without it still work.

## Decisions for you

- **D1. Stop tracking generated files?** Recommended: `git rm -r --cached build_output data/backup`, add `/build_output/` and `/data/backup/` to `.gitignore`, and delete the stray `gitignore` file. This changes repo history going forward, so it's your call.
**Decided 2026-11-11 (done on this branch):**
- **D2:** relabeled "How standard" to **"Campaign occurrence"**. The Definitions panel now says plainly that it measures the share of campaigns studying the assay's *function*, and that sibling formats show the same value. The technology-specific option below remains open for later.
- **D3:** added an **"Include unspecified results"** switch under Inhibition / Activation (default off = only campaigns labeled with the chosen direction, or mixed). The pipeline writes two extra effect views (`inhibition_labeled`, `activation_with_unspecified`); the old `inhibition` / `activation` rows keep their meaning. The title-keyword key for better labels is written out in [plans/EFFECT_LABELING.md](plans/EFFECT_LABELING.md); measured, it changes 4.0% of records.
- **D5:** class-wide fallbacks and small samples now show a **⚠** next to the occurrence. Its tooltip (hover or keyboard focus) says how many campaigns the protein itself has; a Definitions entry explains it.
- **D6:** plan and proposed tightenings and expansions in [plans/REFACTOR_function_classes.md](plans/REFACTOR_function_classes.md). Nothing in the curated files was edited.
- **D8:** `assay_db.json` usage rows are now a compact table (schema 0.2) with the shared source text stored once. `index.html` expands it on load (and still reads the old format). **42,123 rows in 4.9 MB, down from 19,555 rows in 12.6 MB.** `data/usage_stats.csv` was regenerated with the duplicate-study fix and the new effect views. The BindingDB links were computed in memory; nothing in `build_cache/` or `build_output/` was written. Backup: `data/backup/usage_stats_20261111_*.csv`.

**Still open:**
- **D2b. Make campaign occurrence format-specific (later).** Add an `evidence_technology` column to `assays.csv` listing the technology labels each entry uses (for example IC3 → "Automated patch clamp"; IC2 → "Thallium flux"; IC1 → "Membrane potential dye"; P1 → "Fluorogenic / FRET substrate"). I would then change `using_for` so a campaign counts only when both its function and its technology match, and show the function-level number as context. Alternative without a data change: relabel the bar "share of campaigns studying this function".
- **D3. Better direction labels.** Adopt the title-first rule BindingDB already uses (with potentiate, PAM, induce and opener added) for ChEMBL and PubChem too. That needs a relabel (the cached embeddings make it mostly fast). Separately, decide whether "unspecified" should stay in the inhibition view.
- **D4. Class attribution.** Count a campaign toward a class only if its function label belongs to that class (for example a kinase-activity campaign does not count for "autophagy").
- **D5. Small-sample display.** Show the class fallback more loudly (for example "protein has 3 campaigns; showing all kinases"), and hide the share when even the best tier has fewer than 5.
- **D6. Tighten class rules in `function_classes.csv`**: protein_kinase keyword → "Serine/threonine-protein kinase|Tyrosine-protein kinase"; ion_channel: drop GO "channel activity" (keep "ion channel activity").
- **D7. BindingDB literature rows.** Keep them (default), or run `bindingdb_label.py --no-literature` until ChEMBL documents are matched by PMID.
- **D8. Shrink `assay_db.json`** (later): move the repeated `source`/`notes` text to a lookup table. This changes the JSON shape, so `index.html` would change with it.

## How to run

All commands are for PowerShell in `C:\Users\xanox\Documents\GitHub\bioassay-builder`.

**Check the tests** (offline, about 1 second):
```powershell
py -m unittest discover -s tests -v
```

**Your data is already updated** (`data\usage_stats.csv`, `data\condition_stats.csv`, `docs\assay_db.json`, done 2026-11-11). Once, before your next pipeline run, store the BindingDB links in the cache so future rebuilds keep the duplicate-study fix, and refresh the files in `build_output\` (no internet; minutes; `data\usage_stats.csv` is backed up first):
```powershell
py scripts\bindingdb_label.py --rekey                         # adds ChEMBL-document links to existing BindingDB summaries
py scripts\build_database.py --aggregate-only --merge         # recount with one key per study; prunes stale rows
py tools\build_db.py                                          # rebuilds docs\assay_db.json
py scripts\equipment_scan.py                                  # optional: corrected equipment counts
py tools\audit_gaps.py                                        # optional: refresh build_output\audit\
```
The numbers should match what is on the site now. If you skip `--rekey`, the next aggregate quietly brings back about 7,000 duplicate campaigns.
To undo the data change, copy the newest file in `data\backup\` over `data\usage_stats.csv` and rerun `py tools\build_db.py`.

**Look at the site locally:**
```powershell
py -m http.server 8000 --directory docs      # then open http://localhost:8000 ; Ctrl+C to stop
```

**Review and commit the branch** (when you're happy) in GitHub Desktop: branch `audit-and-gaps`. Suggested separate commits: (1) pipeline fixes + tests, (2) website fixes, (3) audit tool + AUDIT/GAPS/README.

The full pipeline order and failure recovery are now in [README.md](README.md#running-the-data-pipeline-windows-powershell).

## Files changed on this branch

- `scripts/build_database.py`: campaign identity helpers, re-keying in aggregate, class merge reuse, stale-row pruning, readout-suggestion fix, `merge_effect` blank fix
- `scripts/bindingdb_label.py`: ChEMBL-document aliases, `--rekey`
- `scripts/plate_formats.py`, `scripts/conditions.py`: use the shared helper (same output)
- `scripts/equipment_scan.py`: whole-campaign scan, protein counts, regex fixes
- `scripts/usage_from_m2.py`: counts campaigns, not assays
- `tools/build_db.py`: faster condition check, size warning
- `docs/index.html`: doctype/viewport, header date, regex escape, null guard, drawer focus and Escape, aria-live, wording
- New: `tests/test_pipeline.py` (21 tests), `tools/audit_gaps.py`, `AUDIT.md`, `GAPS.md`, README pipeline section, `build_output/audit/*.csv`
