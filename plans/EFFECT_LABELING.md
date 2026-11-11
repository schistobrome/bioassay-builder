# Proposal: direction labels from title keywords (inhibition / activation / mixed / unspecified)

Status: **proposal, not implemented**. The website now has an "Include unspecified results" switch as a stopgap (decision D3, 2026-11-11). This document says how to label direction properly, using the approach that already works for BindingDB (`scripts/bindingdb_local.py`, `effect_of`), extended and written out as an explicit key.

## Why change

Today `scripts/m2_label.py` looks for "inhibit" anywhere in the *whole* text first, then "activat/agonist/EC50/stimulat/enhance". Protocol text is full of stray words ("IBMX, a phosphodiesterase inhibitor", "heat-inactivated serum", "agonist-induced"), and the rule misses positive allosteric modulators, potentiators, openers and correctors. BindingDB rows are already labeled from the assay *name* (the first 250 characters), which is much cleaner.

Measured on 57,013 labeled ChEMBL + PubChem records (800 random proteins), the key below changes **2,256 records (4.0%)**:

| current label → key label | records | typical example |
|---|---|---|
| unspecified → activation | 1,037 | "Positive allosteric modulator activity at human mGlu6 receptor … thallium flux" |
| inhibition → activation | 457 | "Agonist activity at mGlu6 … assessed as inhibition of forskolin-stimulated cAMP" (Gi agonist) |
| inhibition → unspecified | 324 | "Biochemical RhoGAP selectivity assay" (no direction in the title; "inhibit" was in the protocol) |
| activation → unspecified | 165 | "hmGluR6_Thallium_Fold_Shift (GIRK)" |
| unspecified → inhibition | 160 | |
| activation → inhibition | 95 | "Inverse agonist activity at human CRTH2 … [35S]GTPγS" |
| inhibition → mixed | 18 | "agonists and antagonists of …" |

The prototype is in the session scratchpad (`title_effect_v3.py`) and is not part of the repo; the key below is complete enough to reimplement it.

## The key

**What is read.** The *head* of the record: at most the first 250 characters of

- ChEMBL: the assay description (the one-line title);
- PubChem: the assay name (`title`);
- BindingDB: the assay description (as today).

The full text (description + protocol) is used only in step 6.

**Step 0. Remove noise phrases** (case-insensitive) before matching. These name reagents or stimuli, not what the compound does:

| Noise phrase (regex idea) | Why it is noise |
|---|---|
| `light activation`, `activated charcoal`, `activated sodium (ortho)vanadate`, `sodium orthovanadate` | reagents |
| `activating mutation(s)`, `activating mutant(s)`, `constitutively activ…`, `activated form`, `activated state`, `activation loop`, `activation segment` | describe the protein |
| `chelation-enhanced`, `enhanced fluorescence`, `enhanced chemiluminescence` | detection chemistry |
| `heat-inactivated` | serum |
| `phosphodiesterase / PDE / protease / phosphatase inhibitor(s) (cocktail)`, `IBMX` | assay additives |
| `<stimulus>-induced / -stimulated / -evoked / -mediated`, where the stimulus is agonist, ligand, forskolin, ATP, glucose, ionomycin, LPS, TNF…, IL-… | the stimulus in an antagonist-mode or cell assay |

**Step 1. Both directions named on purpose → `mixed`.** `agonist(s) and/or antagonist(s)`, `activator(s) and/or inhibitor(s)` (either order, joined by "and", "or" or "/").

**Step 2. Pharmacology role words.** These say what the compound *is*, so they win over read-out verbs. If both lists match, the word that comes **first** wins.

| → inhibition | → activation |
|---|---|
| antagonist, antagonism, antagonistic | agonist, agonism, agonistic (not preceded by "inverse"), partial agonist, full agonist |
| inverse agonist | positive allosteric modulator / modulation |
| negative allosteric modulator / modulation | potentiator |
| blocker | activator |
| inhibitor | opener (channel opener) |
| degrader | corrector (e.g. CFTR) |

Do **not** use the bare abbreviations NAM and PAM: "NAM" is nicotinamide in NAMPT assays.

**Step 3. Gi read-out rule.** `inhibition of forskolin-induced/-stimulated (cAMP)` with no role word → **activation** (Gi-coupled agonist).

**Step 4. Read-out verbs.** Again, the first match wins.

| → inhibition | → activation |
|---|---|
| inhibit, inhibition, inhibitory | activation, activate(s/d/ing) |
| inactivation, inactivate (covalent kinact) | stimulation, stimulate |
| block, blockade, blocking | enhancement, enhance |
| suppress, suppression | induction, induce |
| attenuate, attenuation | potentiation, potentiate |
| deactivate | upregulation |
| degradation, reduction of | stabilization (e.g. pharmacological chaperone) |

**Step 5. Data type, when no word matched.** `IC50` in the head → inhibition; `EC50` → activation. For BindingDB, also the affinity columns, as today: IC50 only → inhibition, EC50 only → activation, Ki/Kd only → unspecified (binding). `AC50` is not a direction.

**Step 6. Binding or no direction.** If the head says binding, displacement, competition, radioligand, affinity, Kd, Ki, thermal shift, target engagement or occupancy → **unspecified**. Otherwise apply steps 2 and 4 to the noise-cleaned *full* text (first 1,500 characters); if still nothing → **unspecified**.

Campaign level (unchanged): a campaign whose records disagree becomes `mixed`; unspecified records do not dilute a direction (`merge_effect`).

## How to implement (when approved)

1. Put the key in **one** function, `direction_of(head, full_text="", affinity_types=())`, in a new small module `scripts/direction.py`, with the tables above as named lists. Use it from `m2_label.label_records` (ChEMBL / PubChem) and from `bindingdb_label.py` (replacing `bindingdb_local.effect_of`).
2. Tests: one synthetic title per table row, plus the regression cases in this document (Gi agonist, inverse agonist, NAMPT "NAM", "agonist-induced" in an antagonist title).
3. Make the direction part of the label cache key (as the term-table hash is), so changing the key relabels direction only. No downloads, and the meaning-model cache is reused.
4. Rebuild: `py scripts\build_database.py --aggregate-only --merge` after relabeling, then `py tools\build_db.py`.
5. Spot-check 100 changed records by hand; aim for at least 90% agreement before switching the website's default.
