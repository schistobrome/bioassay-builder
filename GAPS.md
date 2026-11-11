# Knowledge gaps: where Assay Builder is thin or unreliable

Audit of 2026-11-10/11, based on `data/` and `build_output/` as committed (19,555 usage rows, 5,852 summarized proteins).
Rerun the numbers any time with `py tools/audit_gaps.py` (see "How to rerun" at the end).

Supporting tables in `build_output/audit/`:

| File | What it holds |
|---|---|
| `gaps_by_class.csv` | per function class: playbook entries, kits, proteome size, proteins with no or few campaigns, direction balance |
| `entries_missing_fields.csv` | per playbook entry: missing evidence labels, DOIs, kits, plate formats, usage rows; entries that always show the same share as a sibling |
| `vendor_sourcing_candidates.csv` | target or class, assay type, why, suggested vendors or sources, priority, and whether I spot-checked it |
| `proteome_zero_campaigns.csv` | every reviewed human protein with no campaign in our data, with its class |
| `unclassified_readouts.csv` | readouts the website's instrument filter cannot place |

## The numbers in one place

- **Proteome coverage.** 20,429 reviewed human proteins. 5,761 (28%) have at least one campaign; **14,668 (72%) have none**. Of those, 14,577 were never fetched because they have no ChEMBL single-protein target. The local PubChem bulk files add only 157 such proteins (283 assays), so this is a real public-data desert, not a pipeline omission. 9,836 of the zero-campaign proteins match no function class at all.
- **Thin proteins.** Of 5,852 summarized proteins, 1,043 have 0 campaigns with a function label and 3,116 have 1 to 4. So **71% of proteins fall below the website's 5-campaign threshold**, and the site silently shows class-wide numbers for them.
- **Playbook (78 entries).** No entry has a `sources_doi` (78/78) or a curated `plate_formats` (78/78). 46 are drafts. 65 have no kit. 6 have no `evidence_labels`. 25 have no usage rows on any protein. `citations.csv` is empty.
- **Kits.** 15 kits in total, covering only GPCR (4), HDAC (4), BRD4 (3), kinase/helicase ATPase (ADP-Glo) and 3 unverified services. **14 of the 21 classes have no kit at all.**
- **Activation.** Only GPCR (30% of campaigns), nuclear receptor (48%) and ion channel (12%) have meaningful activation data. 12 classes have less than 3%, and epigenetic writers, demethylases and readers have 6, 4 and 5 activation campaigns in total.
- **Duplicate studies (fixed in code on this branch, not yet in the committed data).** Counts are inflated by copies of the same patent study. When the fix is applied, campaigns drop 2–14% per class (kinases 15,307 → 13,232) and 122 usage shares move by more than 5 points. See AUDIT.md issue 1.

## Top 15 gaps, ranked

Severity means how much it can mislead a scientist choosing an assay today. Sources marked ✔ were spot-checked on the vendor's site (see "What I checked on the web").

### 1. Headline "how standard" share cannot tell sibling assays apart. Severity: high
- **Missing:** a technology link per entry. Usage is computed from the *function* label, so 21 of 72 entries always show the same headline number as a sibling: IC1/IC2/IC3, P1/P2, T1/T2, T3/T4, H2/H3, N1/N2, Ph1/Ph2, PPI1/PPI3, RNA1/RNA2, TOP1/TOP2. On KCNH2 (hERG), membrane-potential dye, thallium flux and patch clamp all show "94%".
- **Why it matters:** the share is meant to answer "is this the standard format?". Today it answers "is this the standard *function*?", which is the same for every format.
- **Fill from:** our own data (technology labels already exist per campaign). Add an `evidence_technology` column to `assays.csv` (decision D2 in AUDIT.md). No outside source needed.

### 2. Kits for ion channels: none recorded. Severity: high
- **Missing:** thallium-flux and membrane-potential kits, stable cell lines, patch-clamp services, for 490 channel proteins and 8,060 campaigns.
- **Why it matters:** hERG (IC5) is a counterscreen in nearly every project, and the cell line is usually the rate-limiting reagent.
- **Fill from:** ✔ Thermo Fisher FluxOR (F10016/F10017), Molecular Devices FLIPR Potassium and Membrane Potential kits, Eurofins Discovery and Charles River cell lines, Sophion/Nanion application notes.

### 3. Degrader (PROTAC) entries rest on no data at all. Severity: high
- **Missing:** `degrader_modality` has 0 campaigns by construction (it is a modality, not a protein class). U3, U4 and U5 have no usage rows and no kits.
- **Why it matters:** degrader projects are common, and the site shows "no campaign data yet" for all three steps.
- **Fill from:** ✔ Promega HiBiT and NanoBRET ternary-complex kits and app notes AN330/AN333; BPS Bioscience CRBN/VHL TR-FRET kits; review articles on degrader assay cascades (PubMed).

### 4. Nuclear receptor coactivator and reporter kits missing. Severity: high
- **Missing:** kits for NR1 (TR-FRET coactivator) and NR2 (reporter) in the class with the best activation data (1,916 activation-relevant campaigns).
- **Why it matters:** NR agonist and antagonist programs both start here, and a scientist needs the exact receptor-specific kit.
- **Fill from:** ✔ Thermo LanthaScreen coactivator kits (FXR A15140, AR A15878, LXRα A15128, CAR A15141, GR A15899); INDIGO Biosciences reporter kits; Thermo GeneBLAzer.

### 5. Histone methyltransferase and demethylase assays: no verified kit. Severity: high
- **Missing:** E1/E2 have only "RapidFire (unverified)". The activation views have 6 and 4 campaigns.
- **Why it matters:** writers and erasers are major target classes, and the universal SAH-detection assays are the obvious HTS starting point.
- **Fill from:** ✔ Promega MTase-Glo (V7601/V7602); BPS chemiluminescent HMT and LSD1/KDM kits; Revvity AlphaLISA/LANCE Ultra histone-mark kits; Cayman.

### 6. Proteases and phosphatases: no kits. Severity: high
- **Missing:** kits for P1 (625 proteases, 6,875 campaigns) and Ph1/Ph2 (196 phosphatases, 2,776 campaigns).
- **Why it matters:** these are routine fluorogenic kits a scientist would buy rather than build.
- **Fill from:** ✔ BPS Bioscience (90+ fluorogenic kits: 3CL, MMPs, caspases, cathepsins); Promega Caspase-Glo/Protease-Glo; Enzo; Thermo EnzChek (DiFMUP); Millipore and R&D malachite-green kits.

### 7. Most proteins fall back to class-wide numbers without a visible warning. Severity: high
- **Missing:** protein-specific data. 71% of summarized proteins have fewer than 5 classifiable campaigns, so the card shows a small "class-wide" pill. The banner only appears when *no* entry has gene-level data.
- **Why it matters:** a class-wide share (for example all 580 kinases) is presented in the same bar as a protein-specific one.
- **Fill from:** mostly a display decision (AUDIT.md D5). More data from patents (SureChEMBL) and literature (ChEMBL documents already used) helps only a little.

### 8. 72% of the proteome has no campaigns; 12,223 proteins match no function class. Severity: high (for coverage)
- **Missing:** both data and a class, so these proteins get only the generic biophysics entries (G1–G7).
- **Why it matters:** a user entering a transcription factor, metabolic enzyme or scaffold gets almost nothing.
- **Fill from:** new function classes (metabolic enzymes such as dehydrogenases, small-molecule kinases, transferases; transcription factors) built from UniProt/Pfam; Assay Guidance Manual chapters; Chemical Probes Portal and Probes & Drugs for which ones have ever been drugged.

### 9. Activation view rests on almost no data outside GPCR, NR and ion channels. Severity: medium-high
- **Missing:** activation campaigns. 12 classes have less than 3% (kinase 4%, epigenetic 0.5–3%, helicase 1%). On top of that, about 2,200 GPCR, NR and ion-channel records whose title says agonist, potentiator or PAM are labeled "unspecified" and counted toward inhibition.
- **Why it matters:** activator projects (kinase activators such as AMPK or GCK, SIRT activators, channel openers) see rough class-level guesses.
- **Fill from:** better labeling of existing data first (AUDIT.md D3). Then the literature on activator assays (AMPK, glucokinase, SIRT1, KCNQ openers) and NCATS qHTS activator AIDs.

### 10. No DOIs, no citations, no review dates. Severity: medium-high
- **Missing:** `sources_doi` is empty for all 78 entries, `citations.csv` is empty, and `reviewed_on` is blank everywhere (the header said "YYYY-MM-DD"; that display is fixed).
- **Why it matters:** a scientist cannot check where a recommendation comes from, and reviewers cannot be credited.
- **Fill from:** Assay Guidance Manual (NCBI Bookshelf) chapter per assay type; NCATS assay-development papers; vendor technical notes; `scripts/find_literature.py` already finds candidates.

### 11. Plate formats are all estimated. Severity: medium
- **Missing:** `plate_formats` is blank for 78/78 entries. The site uses campaign text where at least 5 campaigns state a format, and otherwise guesses from the HTS rating.
- **Why it matters:** 1536-well compatibility decides instrument and automation choices.
- **Fill from:** vendor datasheets (most list validated formats); `usage_stats.csv` plate rows as a first draft for the reviewer.

### 12. Kinases: one kit for 580 proteins; no binding or selectivity service. Severity: medium
- **Missing:** a second biochemical format (TR-FRET/HTRF), binding kits, and the K6 selectivity panel services.
- **Why it matters:** K1 confirmation and K6 profiling are standard purchases.
- **Fill from:** Thermo LanthaScreen and Z'-LYTE; Revvity HTRF KinEASE; Reaction Biology HotSpot; Eurofins DiscoverX KINOMEscan; Promega NanoBRET TE.

### 13. PPI, E3 ligase, DUB and RNA-binding classes: huge, sparsely covered. Severity: medium
- **Missing:** kits and data. PPI has 1,527 proteins (1,060 with no campaigns), E3 427 (331), RNA-binding 1,653 (831, and 695 of the 824 summarized have fewer than 5 usable campaigns), DUB 117. E3 and DUB have a single entry each.
- **Why it matters:** recommendations here are generic and unbacked.
- **Fill from:** Revvity AlphaLISA/HTRF PPI kits, Promega NanoBiT, BPS (PD-1/PD-L1, KRAS-SOS1), LifeSensors and Boston Biochem ubiquitin reagents, Promega DUB-Glo, RBP-screening reviews.

### 14. Readouts the data uses that the playbook does not list. Severity: medium
- **Missing:** in `readout_suggestions.csv` (after this branch's fix removes false flags for generic biophysics), the largest uncovered technologies are Western blot (6,592 campaign-pairs), fluorogenic/FRET substrate (6,060), radioligand binding (4,964), radiometric (3,873), luciferase reporter (2,087), KINOMEscan (2,044) and mobility shift (1,423) across many classes. Part of this is attribution noise: a protein in two classes lends its kinase campaigns to the other class (AUDIT.md issue 9).
- **Why it matters:** these are real formats (for example radiometric kinase, mobility shift) missing from readout lists, so the instrument filter cannot dim correctly.
- **Fill from:** the reviewer. Start from the per-class rows with the largest `n_campaigns` and `function_label` equal to that class's own function.

### 15. Class rules are too broad for kinases and ion channels. Severity: medium
- **Missing:** precise class rules. The bare `kinase` keyword pulls in metabolic kinases (PGK2, MVK, glycerol kinases) and even chemokines (CCL3, CCL8) as "protein_kinase". The GO term `channel activity` pulls aquaporins, claudins and UCP1 into "ion_channel".
- **Why it matters:** those proteins get kinase or ion-channel assay recommendations that do not apply.
- **Fill from:** UniProt keyword list ("Serine/threonine-protein kinase", "Tyrosine-protein kinase") and InterPro family IDs. This is a change to `function_classes.csv` for the reviewer (AUDIT.md D6).

## Smaller items

- 3 playbook readouts the instrument filter cannot place: IC2 "rubidium efflux with atomic absorption", IC5 "binding", E9 "Competition binding panel". Every technology label in the data *is* placed.
- The equipment list in `equipment_mentions.csv` had brand false positives (fixed on this branch; rerun `equipment_scan.py`).
- GPCR orphans and adhesion GPCRs: 480 GPCRs have no campaigns, but 439 of them are olfactory or taste receptors. The rest (ADGR*, FZD3/6, GPR*) are a real but small gap.

## What I checked on the web (2026-11-10)

I used web search result summaries (page text extracted by the search tool). I did not open full datasheets or PDFs, so catalog numbers and formats should still be confirmed before anything goes into `kits.csv`. I added nothing to `kits.csv`.

| Gap | What the vendor offers | Source |
|---|---|---|
| Ion channels (2) | FluxOR thallium-flux potassium channel kit; F10016 (10 plates) and F10017 (100 plates); covers hERG and Kir2.1; optional hERG BacMam reagent | https://thermofisher.com/order/catalog/product/F10016 |
| Nuclear receptors (4) | LanthaScreen TR-FRET coactivator kits per receptor (Tb anti-GST, GST-LBD, fluorescein peptide; 384-well; agonist and antagonist mode). Seen: FXR A15140, AR A15878, LXRα A15128, CAR A15141, GR manual A15899 | https://www.thermofisher.com/order/catalog/product/mx/en/A15878 ; https://www.fishersci.com/shop/products/lanthascreen-tr-fret-fxr-coactivator-assay-kit-rabbit-1/A15140 |
| Degraders (3) | NanoBRET ternary complex (VHL / CRBN) starter kits; HiBiT endogenous tagging for DC50/Dmax; NanoBRET ubiquitination (AN333) | https://at.promega.com/products/protein-detection/protein-degradation-protacs/ ; https://www.promega.com/-/media/files/resources/application-notes/nanobret/an330.pdf |
| HMT (5) | MTase-Glo (V7601/V7602): universal SAH-to-light assay for lysine and arginine HMTs, DNMTs and others; signal half-life over 4 h | https://promega.com/products/epigenetics/methylation-analysis/mtase_glo-methyltransferase-assay |
| Proteases (6) | BPS Bioscience lists 90+ fluorogenic kits (3CL, MMPs, caspases, TMPRSS2, cathepsin L); 96- and 384-well formats | https://bpsbioscience.com/product-types/biochemical-assay-kits-by-format-type/fluorogenic |

## How to rerun

```powershell
cd "C:\Users\xanox\Documents\GitHub\bioassay-builder"
py tools\audit_gaps.py
# whole-proteome numbers need UniProt's table once (about 3 MB):
Invoke-WebRequest -OutFile "$env:TEMP\human_reviewed.tsv.gz" -Uri "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&query=(organism_id:9606)%20AND%20(reviewed:true)&fields=accession,gene_primary,protein_name,length,ec,keyword,go_f,ft_domain,ft_topo_dom,ft_region,ft_motif,ft_repeat,ft_zn_fing,ft_dna_bind"
py tools\audit_gaps.py --uniprot-tsv "$env:TEMP\human_reviewed.tsv.gz"
```
