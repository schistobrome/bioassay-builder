"""
Small offline tests for the pipeline (synthetic data only, no internet, nothing in build_cache/ is touched).

Run from the repo folder:
    py -m unittest discover -s tests -v

Each test builds a tiny fake cache in a temporary folder, points the scripts at it, and checks the result.
"""
import argparse
import csv
import json
import os
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import build_database as bd          # noqa: E402


def write_json(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle)


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def camp(key, f=("Kinase activity",), t=(), e="inhibition", x="cell-free"):
    return {"c": key, "f": list(f), "t": list(t), "e": e, "x": x, "r": ["unspecified"]}


ASSAYS = [
    {"entry_id": "K1", "classes": "protein_kinase", "scope": "cell-free", "evidence_labels": "Kinase activity",
     "readouts": "TR-FRET; luminescence"},
    {"entry_id": "K2", "classes": "protein_kinase", "scope": "cell-based", "evidence_labels": "Kinase activity",
     "readouts": "cellular phosphorylation"},
    {"entry_id": "G1", "classes": "generic", "scope": "cell-free", "evidence_labels": "", "readouts": "DSF"},
]


class FakeCache(unittest.TestCase):
    """Points build_database.py at an empty temporary folder for each test."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="bd_test_")
        self.saved = {k: getattr(bd, k) for k in ("CACHE", "OUTDIR", "SUMMARY", "SUMMARY_BDB", "PLATES_DIR",
                                                   "CONDS_DIR", "DATA")}
        bd.CACHE = os.path.join(self.tmp, "build_cache")
        bd.OUTDIR = os.path.join(self.tmp, "build_output")
        bd.SUMMARY = os.path.join(bd.CACHE, "summary")
        bd.SUMMARY_BDB = os.path.join(bd.CACHE, "summary_bdb")
        bd.PLATES_DIR = os.path.join(bd.CACHE, "plates")
        bd.CONDS_DIR = os.path.join(bd.CACHE, "conditions")
        bd.DATA = os.path.join(self.tmp, "data")
        for d in (bd.SUMMARY, bd.OUTDIR, bd.DATA):
            os.makedirs(d, exist_ok=True)

    def tearDown(self):
        for k, v in self.saved.items():
            setattr(bd, k, v)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def summary(self, acc, gene, campaigns, classes=("protein_kinase",)):
        write_json(os.path.join(bd.SUMMARY, acc + ".json"),
                   {"acc": acc, "gene": gene, "scope": gene, "classes": list(classes),
                    "n_campaigns": len(campaigns), "campaigns": campaigns})

    def run_aggregate(self, **flags):
        args = argparse.Namespace(min_readout=2, suggest_min=1, no_bindingdb=False, **flags)
        return bd.aggregate(args, ASSAYS, [], "2026-01-01")

    def usage(self):
        return read_csv(os.path.join(bd.OUTDIR, "usage_stats.generated.csv"))


# ---------------------------------------------------------------------------------------------- campaign identity

class TestCampaignIdentity(unittest.TestCase):

    def test_chembl_copy_joins_pubchem_original(self):
        records = [
            {"id": "pubchem:100", "source": "PubChem/NCGC", "campaign": "pc:100", "pubchem_aid": "100",
             "function_labels": "Kinase activity", "technology_labels": "", "effect": "inhibition", "context": "cell-free"},
            {"id": "chembl:CHEMBL1", "source": "ChEMBL", "campaign": "doc:D9", "pubchem_aid": "100",
             "function_labels": "", "technology_labels": "TR-FRET", "effect": "unspecified", "context": "cell-free"},
            {"id": "chembl:CHEMBL2", "source": "ChEMBL", "campaign": "doc:D5", "pubchem_aid": "",
             "function_labels": "Kinase activity", "technology_labels": "", "effect": "activation", "context": "cell-based"},
        ]
        campaigns, merged = bd.dedupe_campaigns(records)
        self.assertEqual(merged, 1)
        by_key = {c["c"]: c for c in campaigns}
        self.assertEqual(sorted(by_key), ["doc:D5", "pc:100"])
        self.assertEqual(by_key["pc:100"]["t"], ["TR-FRET"])           # labels from both copies are pooled
        self.assertEqual(by_key["pc:100"]["e"], "inhibition")           # unspecified does not dilute a direction
        self.assertEqual(by_key["doc:D5"]["x"], "cell-based")

    def test_merge_effect(self):
        self.assertEqual(bd.merge_effect(["inhibition", "unspecified"]), "inhibition")
        self.assertEqual(bd.merge_effect(["inhibition", "activation"]), "mixed")
        self.assertEqual(bd.merge_effect(["unspecified", ""]), "unspecified")    # a blank is not a direction
        self.assertEqual(bd.merge_effect([]), "unspecified")

    def test_effect_subset(self):
        cs = [camp("a", e="inhibition"), camp("b", e="activation"), camp("c", e="mixed"), camp("d", e="unspecified")]
        self.assertEqual([c["c"] for c in bd.effect_subset(cs, "inhibition")], ["a", "c", "d"])
        self.assertEqual([c["c"] for c in bd.effect_subset(cs, "activation")], ["b", "c"])
        self.assertEqual([c["c"] for c in bd.effect_subset(cs, "inhibition_labeled")], ["a", "c"])
        self.assertEqual([c["c"] for c in bd.effect_subset(cs, "activation_with_unspecified")], ["b", "c", "d"])

    def test_canonical_patent_keys(self):
        self.assertEqual(bd.canonical_campaign("patent:US9023767B2"), "patent:US9023767")
        self.assertEqual(bd.canonical_campaign("patent:US9023767"), "patent:US9023767")
        self.assertEqual(bd.canonical_campaign("patent:US20240368106A1"), "patent:US20240368106")
        self.assertEqual(bd.canonical_campaign("patent:WO2012123456"), "patent:WO2012123456")
        self.assertEqual(bd.canonical_campaign("doc:CHEMBL1"), "doc:CHEMBL1")
        self.assertEqual(bd.canonical_campaign("doc:CHEMBL1", {"doc:CHEMBL1": "patent:US1234567"}), "patent:US1234567")


# ---------------------------------------------------------------------------------------------- text extraction

class TestTextExtraction(unittest.TestCase):

    def test_plates(self):
        self.assertEqual(bd.plates_in("Assays were run in 384-well plates."), ["384"])
        self.assertEqual(bd.plates_in("1536 well format; compounds from a 384-well source plate"), ["1536"])
        self.assertEqual(bd.plates_in("cells seeded in 96 well microtiter plates and 6-well dishes"), ["96", "low"])
        self.assertEqual(bd.plates_in("IC50 of 384 nM"), [])
        self.assertEqual(bd.plates_in("compound plates (384-well) were stamped into 1536-well assay plates"), ["1536"])

    def test_conditions(self):
        text = ("Reactions in 50 mM HEPES pH 7.5, 10 mM MgCl2, 1 mM DTT, 0.01% Tween-20, 0.1% BSA, 10 uM ATP, "
                "1% DMSO final, incubated for 60 min at room temperature; final volume 10 uL.")
        c = bd.conditions_in(text)
        self.assertEqual(c["ph"], 7.5)
        self.assertEqual(c["mgcl2_mm"], 10.0)
        self.assertEqual(c["atp_um"], 10.0)
        self.assertEqual(c["dmso_pct"], 1.0)
        self.assertEqual(c["bsa_pct"], 0.1)
        self.assertEqual(c["incubation_min"], 60.0)
        self.assertEqual(c["volume_ul"], 10.0)
        self.assertEqual(c["temperature"], ["room temperature"])
        self.assertEqual(c["buffer"], ["HEPES"])
        self.assertEqual(c["detergent"], ["Tween-20"])
        self.assertEqual(c["reducing_agent"], ["DTT"])

    def test_conditions_rejects_implausible(self):
        c = bd.conditions_in("Compounds dissolved in 100% DMSO. ATP 1 mM. Incubated 2 h at 37 C. pH 14")
        self.assertNotIn("dmso_pct", c)                      # 100 % is the stock, not the assay
        self.assertEqual(c["atp_um"], 1000.0)
        self.assertEqual(c["incubation_min"], 120.0)
        self.assertEqual(c["temperature"], ["37 C"])
        self.assertNotIn("ph", c)


# ---------------------------------------------------------------------------------------------- aggregate

class TestAggregate(FakeCache):

    def test_gene_rows_and_cell_scope(self):
        self.summary("P1", "KIN1", [
            camp("doc:D1", t=["TR-FRET"]),
            camp("pc:100", t=["TR-FRET"], e="unspecified"),
            camp("pc:101", t=["Luminescence"], e="activation"),
            camp("doc:D2", e="inhibition", x="cell-based"),
        ])
        res = self.run_aggregate()
        rows = self.usage()
        top = {(r["entry_id"], r["effect"]): r for r in rows if r["scope_type"] == "gene" and not r["readout"]}
        k1 = top[("K1", "inhibition")]
        self.assertEqual((k1["n_campaigns_total"], k1["n_campaigns_classifiable"], k1["n_campaigns_using"]),
                         ("3", "3", "2"))                   # the cell-based campaign is not counted for a cell-free entry
        k2 = top[("K2", "inhibition")]
        self.assertEqual(k2["n_campaigns_using"], "1")
        self.assertEqual(top[("K1", "activation")]["n_campaigns_using"], "1")
        self.assertEqual(top[("K1", "inhibition_labeled")]["n_campaigns_total"], "2")            # pc:100 (unspecified) left out
        self.assertEqual(top[("K1", "activation_with_unspecified")]["n_campaigns_total"], "2")
        ids = [r["stat_id"] for r in rows]
        self.assertEqual(len(ids), len(set(ids)))                                                # no stat_id collisions
        tr = [r for r in rows if r["entry_id"] == "K1" and r["effect"] == "inhibition" and r["readout"] == "TR-FRET"
              and r["scope_type"] == "gene"]
        self.assertEqual((tr[0]["n_campaigns_classifiable"], tr[0]["n_campaigns_using"]), ("2", "2"))
        self.assertEqual(res["entries_without_evidence_labels"], ["G1"])
        cls = [r for r in rows if r["scope_type"] == "function_class" and r["entry_id"] == "K1"
               and r["effect"] == "inhibition" and not r["readout"]]
        self.assertEqual(cls[0]["n_campaigns_using"], "2")

    def test_generic_entries_cover_readout_suggestions(self):
        self.summary("P1", "KIN1", [camp("doc:D1", t=["Thermal shift (DSF)"]), camp("doc:D2", t=["Radiometric"])])
        self.run_aggregate()
        sug = {r["technology"]: r for r in read_csv(os.path.join(bd.OUTDIR, "readout_suggestions.csv"))}
        self.assertEqual(sug["Thermal shift (DSF)"]["suggestion"], "covered")           # G1 (generic) lists DSF
        self.assertEqual(sug["Radiometric"]["suggestion"], "NOT in any playbook readout list for this class")

    def test_class_rows_count_a_shared_campaign_once(self):
        self.summary("P1", "KIN1", [camp("pc:100"), camp("doc:D1")])
        self.summary("P2", "KIN2", [camp("pc:100"), camp("doc:D3")])
        self.run_aggregate()
        cls = [r for r in self.usage() if r["scope_type"] == "function_class" and r["entry_id"] == "K1"
               and r["effect"] == "inhibition" and not r["readout"]]
        self.assertEqual(cls[0]["n_campaigns_total"], "3")

    def test_patent_copies_are_one_campaign(self):
        """A PubChem patent deposit (key with kind code) and the BindingDB row for the same patent are one campaign."""
        self.summary("P1", "KIN1", [camp("patent:US9023767B2"), camp("doc:D1")])
        write_json(os.path.join(bd.SUMMARY_BDB, "P1.json"),
                   {"acc": "P1", "campaigns": [camp("patent:US9023767", t=["TR-FRET"]), camp("pmid:55")]})
        self.run_aggregate()
        top = [r for r in self.usage() if r["scope_type"] == "gene" and r["entry_id"] == "K1"
               and r["effect"] == "inhibition" and not r["readout"]]
        self.assertEqual(top[0]["n_campaigns_total"], "3")       # not 4

    def test_bindingdb_alias_joins_chembl_document(self):
        """BindingDB tells us which ChEMBL document is the same patent; the three copies collapse to one."""
        self.summary("P1", "KIN1", [camp("patent:US1111111B2"), camp("doc:CHEMBL7"), camp("doc:D1")])
        write_json(os.path.join(bd.SUMMARY_BDB, "P1.json"),
                   {"acc": "P1", "campaigns": [camp("patent:US1111111")],
                    "aliases": {"doc:CHEMBL7": "patent:US1111111"}})
        self.run_aggregate()
        top = [r for r in self.usage() if r["scope_type"] == "gene" and r["entry_id"] == "K1"
               and r["effect"] == "inhibition" and not r["readout"]]
        self.assertEqual(top[0]["n_campaigns_total"], "2")
        effect = read_csv(os.path.join(bd.OUTDIR, "effect_by_gene.csv"))
        self.assertEqual(effect[0]["n_campaigns"], "2")

    def test_no_bindingdb_flag_ignores_aliases(self):
        self.summary("P1", "KIN1", [camp("doc:CHEMBL7"), camp("doc:D1")])
        write_json(os.path.join(bd.SUMMARY_BDB, "P1.json"),
                   {"acc": "P1", "campaigns": [camp("patent:US1111111")], "aliases": {"doc:CHEMBL7": "patent:US1111111"}})
        args = argparse.Namespace(min_readout=2, suggest_min=1, no_bindingdb=True)
        bd.aggregate(args, ASSAYS, [], "2026-01-01")
        top = [r for r in self.usage() if r["scope_type"] == "gene" and r["entry_id"] == "K1"
               and r["effect"] == "inhibition" and not r["readout"]]
        self.assertEqual(top[0]["n_campaigns_total"], "2")

    def test_plates_attach_to_original_key(self):
        self.summary("P1", "KIN1", [camp("patent:US1111111B2"), camp("doc:D1"), camp("doc:D2"), camp("doc:D3"),
                                    camp("doc:D4")])
        write_json(os.path.join(bd.PLATES_DIR, "P1.json"),
                   {"campaigns": {"patent:US1111111B2": ["384"], "doc:D1": ["384"], "doc:D2": ["1536"]}})
        self.run_aggregate()
        plates = {r["readout"]: r for r in self.usage() if r["scope_type"] == "gene" and r["entry_id"] == "K1"
                  and r["effect"] == "inhibition" and r["readout"].startswith("plate:")}
        self.assertEqual(plates["plate:384"]["n_campaigns_using"], "2")
        self.assertEqual(plates["plate:384"]["n_campaigns_classifiable"], "3")


class TestBindingDBAliases(unittest.TestCase):

    def test_one_to_one_text_match(self):
        import bindingdb_label as bdl
        long_a = "Kinase Assay: The inhibitory activity of the compounds against JAK1 was measured by HTRF in 384-well plates"
        long_b = "Cell assay: phosphorylation of STAT5 in TF-1 cells stimulated with GM-CSF was measured by AlphaLISA SureFire"
        generic = "Inhibition assay"
        chembl = [
            {"id": "chembl:A1", "source": "ChEMBL", "campaign": "doc:CHEMBL7", "title": long_a},
            {"id": "chembl:A2", "source": "ChEMBL", "campaign": "doc:CHEMBL8", "title": long_b},
            {"id": "chembl:A3", "source": "ChEMBL", "campaign": "doc:CHEMBL9", "title": long_b},   # same text, 2 docs
            {"id": "chembl:A4", "source": "ChEMBL", "campaign": "doc:CHEMBL10", "title": generic},
        ]
        bdb = [
            {"id": "bdb:1", "campaign": "patent:US1111111", "title": long_a[:150]},
            {"id": "bdb:2", "campaign": "patent:US2222222", "title": long_b[:150]},   # ambiguous: two ChEMBL docs
            {"id": "bdb:3", "campaign": "patent:US3333333", "title": generic},         # too short to trust
        ]
        self.assertEqual(bdl.chembl_aliases(bdb, chembl), {"doc:CHEMBL7": "patent:US1111111"})

    def test_rekey_runs_without_folder_and_labeling_still_needs_it(self):
        import bindingdb_label as bdl
        from unittest import mock
        with mock.patch.object(bdl, "rekey_only", return_value="rekeyed") as rk:
            self.assertEqual(bdl.main(["--rekey"]), "rekeyed")
            rk.assert_called_once()
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            bdl.main([])

    def test_chembl_document_linked_to_two_patents_is_not_merged(self):
        import bindingdb_label as bdl
        text = "Biochemical assay: inhibition of recombinant human PDE4B catalytic domain measured by IMAP FP readout"
        chembl = [{"id": "chembl:A1", "source": "ChEMBL", "campaign": "doc:CHEMBL7", "title": text}]
        bdb = [{"id": "bdb:1", "campaign": "patent:US1", "title": text}, {"id": "bdb:2", "campaign": "patent:US2", "title": text}]
        self.assertEqual(bdl.chembl_aliases(bdb, chembl), {})


class TestEquipmentScan(unittest.TestCase):

    def test_whole_campaign_is_read_and_proteins_counted(self):
        import equipment_scan as es
        p1 = [
            {"id": "chembl:C1", "source": "ChEMBL", "campaign": "doc:D1", "pubchem_aid": "100", "text": "Inhibition of X"},
            {"id": "pubchem:100", "source": "PubChem/NCGC", "campaign": "pc:100", "pubchem_aid": "100",
             "text": "Reagents were dispensed with a Multidrop Combi and read on an EnVision plate reader."},
        ]
        p2 = [{"id": "pubchem:100", "source": "PubChem/NCGC", "campaign": "pc:100", "pubchem_aid": "100",
               "text": "Reagents were dispensed with a Multidrop Combi."}]
        hits, cands, n_camp, n_prot, n_text = es.scan([("P1", p1), ("P2", p2)])
        md = hits[("bulk reagent dispenser", "Multidrop / Multidrop Combi")]
        self.assertEqual(md["camps"], {"pc:100"})                  # one campaign, although two proteins and two records
        self.assertEqual(md["prots"], {"P1", "P2"})
        self.assertIn(("plate readers / detectors", "EnVision (PerkinElmer)"), hits)
        self.assertEqual((n_camp, n_prot), (1, 2))

    def test_false_positives_removed(self):
        import equipment_scan as es
        text = ("Geminin degradation in mosquito C6/36 cells; cells were harvested; omega-3 fatty acids; "
                "drug synergy; antagonist at the M5 muscarinic receptor; a Hamilton syringe.")
        hits, *_ = es.scan([("P1", [{"id": "x", "source": "PubChem/A", "campaign": "pc:1", "text": text}])])
        self.assertEqual(sorted(name for _, name in hits), [])


class TestUsageFromM2(unittest.TestCase):

    def test_assays_of_one_paper_count_once(self):
        import usage_from_m2 as um
        records = [{"campaign": "doc:D1", "function_labels": "Kinase activity", "technology_labels": "TR-FRET"},
                   {"campaign": "doc:D1", "function_labels": "Kinase activity", "technology_labels": "Luminescence"},
                   {"campaign": "pc:5", "function_labels": "", "technology_labels": ""}]
        camps = um.campaigns_of(records)
        self.assertEqual(len(camps), 2)
        self.assertEqual(camps[0]["technology_labels"], "Luminescence; TR-FRET")


# ---------------------------------------------------------------------------------------------- merge into data/

class TestMerge(FakeCache):

    COLS = bd.USAGE_COLUMNS

    def row(self, entry, scope, effect, version, readout="", using="1"):
        return {"stat_id": "%s-%s-%s" % (scope, entry, effect), "entry_id": entry, "scope_type": "gene",
                "scope_id": scope, "stage": "any", "effect": effect, "readout": readout, "n_campaigns_total": "5",
                "n_campaigns_classifiable": "5", "n_campaigns_using": using, "method_version": version}

    def test_merge_rules(self):
        self.summary("P1", "KIN1", [camp("doc:D1")])
        self.summary("P2", "KIN2", [camp("doc:D2")])
        existing = [
            self.row("K1", "KIN1", "inhibition", "bd-20260101", using="1"),     # auto, regenerated: replaced
            self.row("K2", "KIN1", "inhibition", "bd-20260101"),                # auto, not produced again: pruned
            self.row("K1", "KIN1", "activation", "hand"),                       # yours: kept even if not regenerated
            self.row("K1", "KIN2", "any", "bd-20260101"),                       # auto, protein summarized, no rows now: pruned
            self.row("K1", "OTHER", "any", "bd-20260101"),                      # auto, protein not summarized: kept
        ]
        bd.write_csv(os.path.join(bd.DATA, "usage_stats.csv"), existing, self.COLS)
        gen = [self.row("K1", "KIN1", "inhibition", "bd-20260202", using="4")]
        bd.write_csv(os.path.join(bd.OUTDIR, "usage_stats.generated.csv"), gen, self.COLS)
        added, replaced, kept, conflicts, pruned = bd.merge_into_data(False)
        final = read_csv(os.path.join(bd.DATA, "usage_stats.csv"))
        keys = sorted((r["entry_id"], r["scope_id"], r["effect"], r["method_version"]) for r in final)
        self.assertEqual(keys, [("K1", "KIN1", "activation", "hand"), ("K1", "KIN1", "inhibition", "bd-20260202"),
                                ("K1", "OTHER", "any", "bd-20260101")])
        self.assertEqual((added, replaced, kept, pruned), (0, 1, 0, 2))
        self.assertTrue(os.listdir(os.path.join(bd.DATA, "backup")))         # a backup was written

    def test_hand_row_wins_over_generated(self):
        bd.write_csv(os.path.join(bd.DATA, "usage_stats.csv"), [self.row("K1", "KIN1", "inhibition", "hand")], self.COLS)
        bd.write_csv(os.path.join(bd.OUTDIR, "usage_stats.generated.csv"),
                     [self.row("K1", "KIN1", "inhibition", "bd-20260202", using="4")], self.COLS)
        added, replaced, kept, conflicts, pruned = bd.merge_into_data(False)
        final = read_csv(os.path.join(bd.DATA, "usage_stats.csv"))
        self.assertEqual([r["method_version"] for r in final], ["hand"])
        self.assertEqual(kept, 1)


# --------------------------------------------------------------------------- labeling fixes (2026-11-11)

class TestMembraneContext(unittest.TestCase):
    """Membrane preparations and GTPgammaS binding are cell-free, even when the receptor was expressed in cells."""

    def ctx(self, text, bao="", cell_type=""):
        return bd.m2.context_of(text, bao, cell_type)

    def test_membrane_assays_are_cell_free(self):
        self.assertEqual(self.ctx("Displacement of [3H]NMS from human M2 receptor expressed in HEK293T cell membranes",
                                  cell_type="HEK293T"), "cell-free")
        self.assertEqual(self.ctx("Agonist activity at human S1P2 receptor expressed in CHO cells assessed as induction of "
                                  "[35S]GTPgammaS binding", cell_type="CHO"), "cell-free")
        self.assertEqual(self.ctx("Agonist-stimulated GTPγS binding", cell_type="CHO"), "cell-free")
        self.assertEqual(self.ctx("Binding affinity to human CB1 receptor", bao="cell membrane format", cell_type="CHO"),
                         "cell-free")

    def test_cell_assays_stay_cell_based(self):
        self.assertEqual(self.ctx("Inhibition of hERG measured by membrane potential dye in HEK293 cells"), "cell-based")
        self.assertEqual(self.ctx("Displacement of [3H]ligand from receptor in intact CHO cells membranes"), "cell-based")
        self.assertEqual(self.ctx("Agonist activity at receptor expressed in CHO cells assessed as cAMP accumulation",
                                  cell_type="CHO"), "cell-based")
        self.assertEqual(self.ctx("Kinase activity of enzyme expressed in Sf9 cells"), "cell-free")       # unchanged rule
        # ChEMBL's loose "expressed in ... cell membrane" on ordinary cell assays (real examples)
        self.assertEqual(self.ctx("Activation of human TRPV4 expressed in HEK293T cell membrane at 30 uM in presence of GSK",
                                  cell_type="HEK293T"), "cell-based")
        self.assertEqual(self.ctx("Antagonist activity at human mGlu1 expressed in AV12 cell membrane assessed as reduction "
                                  "in glutamate-induced calcium influx", cell_type="AV12"), "cell-based")
        self.assertEqual(self.ctx("Antagonist activity at human S1P2 receptor expressed in CHO cell membrane by FLIPR assay",
                                  cell_type="CHO"), "cell-based")
        self.assertEqual(self.ctx("Corrector activity at CFTR expressed in CFBE41o- cells assessed as restoration of "
                                  "plasma membrane density", bao="cell membrane format", cell_type="CFBE41o-"), "cell-based")
        self.assertEqual(self.ctx("Cell membrane integrity measured by LDH release in HepG2 cells"), "cell-based")
        # GTPgammaS binding is a membrane assay even when the text says "uptake"
        self.assertEqual(self.ctx("Activity at S1P2 receptor expressed in CHO cells measured as S1P-induced "
                                  "[35S]GTP-gamma-S uptake", cell_type="CHO"), "cell-free")

    def test_rules_version_is_part_of_the_label_cache_key(self):
        self.assertTrue(bd.label_hash().endswith(bd.m2.RULES_VERSION))


class TestTermTableFixes(unittest.TestCase):
    """Reads scripts/term_table.csv (read-only) and checks the 2026-11-11 synonym fixes."""

    @classmethod
    def setUpClass(cls):
        cls.terms = bd.m2.load_terms()

    def hits(self, axis, text):
        return set(bd.m2.rule_match(self.terms, axis, bd.m2.norm(text)))

    def test_tritiated_displacement_is_radioligand_binding(self):
        text = "Displacement of [3H]-(S)-ligand from human OX1R expressed in CHOK1 cells"
        self.assertIn("Radioligand binding", self.hits("technology", text))
        self.assertIn("Receptor radioligand binding", self.hits("function", text))
        self.assertIn("Radioligand binding", self.hits("technology", "Displacement of [125I]-CXCL12 from CXCR4"))
        self.assertNotIn("Radioligand binding", self.hits("technology", "[3H]thymidine incorporation in cells"))

    def test_ip_one_has_a_function_label(self):
        self.assertIn("Inositol phosphate accumulation", self.hits("function", "Agonist activity by IP-One HTRF assay"))
        self.assertIn("Inositol phosphate accumulation", self.hits("function", "assessed as IP1 accumulation"))
        self.assertIn("Inositol phosphate accumulation",
                      self.hits("function", "assessed as stimulation of [3H]inositol phosphates accumulation"))
        # enzyme assays on inositol phosphates must not get this cell-only label
        self.assertNotIn("Inositol phosphate accumulation",
                         self.hits("function", "Inhibition of SHIP1 inositol phosphate 5-phosphatase activity"))
        self.assertNotIn("Inositol phosphate accumulation",
                         self.hits("function", "Inhibition of inositol trisphosphate 3-kinase ITPKA"))

    def test_discoverx_is_not_kinomescan(self):
        t = self.hits("technology", "PathHunter beta-arrestin assay (DiscoverX)")
        self.assertNotIn("KINOMEscan competition binding", t)
        self.assertIn("Enzyme fragment complementation", t)
        self.assertIn("Enzyme fragment complementation", self.hits("technology", "cAMP measured by HitHunter assay"))
        self.assertIn("KINOMEscan competition binding", self.hits("technology", "DiscoverX KINOMEscan panel"))

    def test_erk_phosphorylation_is_cellular_signaling(self):
        self.assertIn("Cellular signaling", self.hits("function", "assessed as ERK1/2 phosphorylation by AlphaScreen SureFire"))


class TestEcImpliedFunction(unittest.TestCase):
    """2026-11-12: when no UniProt domain names a function term, the EC number may (one clear term, or nothing)."""
    NAMES = {"Kinase activity", "Protease activity", "Exonuclease / nuclease activity", "Helicase unwinding",
             "Cytochrome P450 activity", "Phosphodiesterase activity", "Enzyme activity",
             "DNA / RNA methyltransferase activity", "Deubiquitinase activity"}

    def test_ec_mapping(self):
        f = bd.m2.ec_implied
        self.assertEqual(f(["1.14.14.1", "1.14.14.55"], self.NAMES), "Cytochrome P450 activity")     # CYP3A4
        self.assertEqual(f(["3.1.4.35"], self.NAMES), "Phosphodiesterase activity")                    # PDE5A
        self.assertEqual(f(["4.2.1.1"], self.NAMES), "Enzyme activity")                                # carbonic anhydrase
        self.assertEqual(f(["3.4.19.12"], self.NAMES), "Deubiquitinase activity")                      # USP7, not protease
        self.assertIsNone(f(["3.6.5.2"], self.NAMES))                       # KRAS: no activity is implied
        self.assertIsNone(f(["7.6.2.2"], self.NAMES))                       # ABC transporter
        self.assertIsNone(f(["2.7.11.1", "3.4.21.5"], self.NAMES))          # two different activities: nothing
        self.assertIsNone(f([], self.NAMES))
        self.assertEqual(f(["2.1.1.37"], self.NAMES), "DNA / RNA methyltransferase activity")        # DNMT1
        self.assertEqual(f(["2.1.1.370"], self.NAMES | {"Histone methyltransferase activity"}),
                         "Histone methyltransferase activity")             # exact match: NSD3 is an HMT, not a DNMT
        self.assertEqual(f(["2.7.1.153"], self.NAMES), "Kinase activity")  # PI3K: lipid kinase
        self.assertIsNone(f(["1.14.14.1"], {"Kinase activity"}))             # term not in the table: nothing

    def test_domain_still_wins_and_ec_only_fills_gaps(self):
        terms = bd.m2.load_terms()
        tmp = tempfile.mkdtemp()
        try:
            def entry(acc, domains, ecs):
                write_json(os.path.join(tmp, "uniprot_%s.json" % acc), {
                    "features": [{"type": "Domain", "description": d} for d in domains],
                    "proteinDescription": {"recommendedName": {"ecNumbers": [{"value": e} for e in ecs]}}})
            entry("K", ["Protein kinase"], ["2.7.11.1"])
            entry("W", ["Helicase ATP-binding", "3'-5' exonuclease"], ["3.6.4.12", "3.1.11.1"])
            entry("C", [], ["4.2.1.1"])
            self.assertEqual(bd.m2.implied_function(terms, tmp, "K")[0], "Kinase activity")
            self.assertIsNone(bd.m2.implied_function(terms, tmp, "W")[0])       # WRN-like: several domains, no guess
            name, detail = bd.m2.implied_function(terms, tmp, "C")
            if "Enzyme activity" in {t.name for t in terms}:
                self.assertEqual((name, detail), ("Enzyme activity", "EC 4.2.1.1"))
            else:
                self.assertIsNone(name)
        finally:
            shutil.rmtree(tmp)


class TestClassRules(FakeCache):
    """Class rules (2026-11-11): needs_domain means domain AND another match; classes are recomputed at aggregate."""
    RECEPTOR = {"class_id": "receptor_ectodomain", "label": "Receptor", "ec_prefixes": "", "keywords": "receptor",
                "go_terms": "receptor activity", "domain_patterns": "extracellular", "needs_domain": "yes",
                "modality": "", "excludes": "", "notes": ""}
    KINASE = {"class_id": "protein_kinase", "label": "Kinase", "ec_prefixes": "2.7.11", "keywords": "",
              "go_terms": "", "domain_patterns": "", "needs_domain": "no", "modality": "", "excludes": "", "notes": ""}

    @staticmethod
    def info(keywords=(), go=(), ecs=(), extracellular=True):
        feats = [{"type": "Topological domain", "description": "Extracellular", "start": 1, "end": 100}] if extracellular else []
        return {"ecs": list(ecs), "keywords": list(keywords), "go": list(go), "features": feats}

    def test_extracellular_region_alone_is_not_a_receptor(self):
        cls = [self.RECEPTOR]
        self.assertNotIn("receptor_ectodomain", bd.m4.classify(self.info(), cls))                     # e.g. an aquaporin
        self.assertIn("receptor_ectodomain", bd.m4.classify(self.info(keywords=["Receptor"]), cls))  # real receptor
        self.assertNotIn("receptor_ectodomain", bd.m4.classify(self.info(keywords=["Receptor"], extracellular=False), cls))

    def test_aggregate_recomputes_classes_from_cached_uniprot(self):
        bd_m1 = bd.m1.DATA
        bd.m1.DATA = os.path.join(bd.CACHE, "m1_data")
        try:
            write_json(os.path.join(bd.m1.DATA, "uniprot_P1.json"),
                       {"proteinDescription": {"recommendedName": {"ecNumbers": [{"value": "2.7.11.1"}]}}})
            # the stored summary says "ion_channel", but the cached UniProt entry makes it a kinase under today's rules
            self.summary("P1", "KIN1", [camp("doc:D1", t=["TR-FRET"])], classes=("ion_channel",))
            args = argparse.Namespace(min_readout=2, suggest_min=1, no_bindingdb=False)
            res = bd.aggregate(args, ASSAYS, [self.KINASE], "2026-01-01")
            self.assertEqual(res["reclassed"], 1)
            rows = self.usage()
            self.assertTrue(any(r["entry_id"] == "K1" and r["scope_type"] == "gene" for r in rows))
        finally:
            bd.m1.DATA = bd_m1


class TestCellFreeFunctions(FakeCache):
    """A campaign with a membrane binding assay AND a cell assay counts for cell-free entries too."""

    def test_dedupe_keeps_cell_free_functions_and_technologies(self):
        recs = [{"id": "a1", "campaign": "doc:D1", "context": "cell-free", "function_labels": "Kinase activity",
                 "technology_labels": "TR-FRET", "effect": "inhibition", "role": ""},
                {"id": "a2", "campaign": "doc:D1", "context": "cell-based", "function_labels": "Cell proliferation / viability",
                 "technology_labels": "Luminescence", "effect": "inhibition", "role": ""}]
        (c,), _ = bd.dedupe_campaigns(recs)
        self.assertEqual(c["x"], "cell-based")                         # unchanged: any cell record makes it cell-based
        self.assertEqual(c["ff"], ["Kinase activity"])
        self.assertEqual(c["tf"], ["TR-FRET"])

    def test_old_summaries_without_ff_behave_as_before(self):
        self.assertEqual(bd.cell_free_functions(camp("x", x="cell-based")), set())
        self.assertEqual(bd.cell_free_functions(camp("x", x="cell-free")), {"Kinase activity"})
        merged = bd.merge_campaign_lists([camp("pc:1", x="cell-based")], [camp("pc:1", x="cell-free")])
        self.assertEqual(merged[0]["x"], "cell-based")
        self.assertEqual(merged[0]["ff"], ["Kinase activity"])         # the cell-free copy keeps its function

    def test_mixed_campaign_counts_for_cell_free_entry_with_cell_free_technologies_only(self):
        mixed = dict(camp("doc:D2", f=["Kinase activity", "Cell proliferation / viability"], t=["TR-FRET", "Luminescence"],
                          x="cell-based"), ff=["Kinase activity"], tf=["TR-FRET"])
        self.summary("P1", "KIN1", [camp("doc:D1", t=["TR-FRET"]), mixed])
        self.run_aggregate()
        rows = self.usage()
        top = {(r["entry_id"], r["effect"]): r for r in rows if r["scope_type"] == "gene" and not r["readout"]}
        self.assertEqual(top[("K1", "inhibition")]["n_campaigns_using"], "2")    # cell-free entry: both campaigns
        self.assertEqual(top[("K2", "inhibition")]["n_campaigns_using"], "1")    # cell-based entry: only the mixed one
        k1_reads = {r["readout"] for r in rows if r["entry_id"] == "K1" and r["scope_type"] == "gene" and r["readout"]
                    and not r["readout"].startswith("plate:")}
        self.assertEqual(k1_reads, {"TR-FRET"})                                  # no cell-assay technology leaks in


def rec(rid, campaign, fn, tech, context="cell-free"):
    return {"id": rid, "source": "ChEMBL", "campaign": campaign, "pubchem_aid": "", "context": context,
            "function_labels": fn, "technology_labels": tech, "effect": "inhibition", "role": ""}


class TestReadoutPairing(FakeCache):
    """A readout counts for an entry only if it was on a record that measured the entry's function."""

    def readouts(self, gene="KIN1", entry="K1"):
        return {r["readout"]: r for r in self.usage() if r["entry_id"] == entry and r["scope_type"] == "gene"
                and r["scope_id"] == gene and r["effect"] == "inhibition" and r["readout"]
                and not r["readout"].startswith("plate:")}

    def test_dedupe_pairs_readouts_with_functions(self):
        (c,), _ = bd.dedupe_campaigns([rec("a1", "doc:D1", "Kinase activity", "TR-FRET / HTRF"),
                                       rec("a2", "doc:D1", "Ligand binding affinity", "KINOMEscan competition binding"),
                                       rec("a3", "doc:D1", "Kinase activity", "Western blot", "cell-based")])
        self.assertEqual(c["tpf"], {"Kinase activity": ["TR-FRET / HTRF"],
                                    "Ligand binding affinity": ["KINOMEscan competition binding"]})
        self.assertEqual(c["tpc"], {"Kinase activity": ["Western blot"]})

    def test_binding_readout_in_same_paper_does_not_count_for_activity_entry(self):
        (c,), _ = bd.dedupe_campaigns([rec("a1", "doc:D1", "Kinase activity", "TR-FRET / HTRF"),
                                       rec("a2", "doc:D1", "Ligand binding affinity", "KINOMEscan competition binding")])
        self.summary("P1", "KIN1", [c])
        self.run_aggregate()
        reads = self.readouts()
        self.assertEqual(set(reads), {"TR-FRET / HTRF"})
        self.assertEqual(reads["TR-FRET / HTRF"]["n_campaigns_classifiable"], "1")

    def test_panels_left_out_of_readout_shares_only(self):
        panel = lambda: bd.dedupe_campaigns([rec("p", "doc:PANEL", "Kinase activity", "Radiometric")])[0][0]
        own = bd.dedupe_campaigns([rec("a", "doc:D1", "Kinase activity", "TR-FRET / HTRF")])[0][0]
        self.summary("P1", "KIN1", [own, panel()])
        self.summary("P2", "KIN2", [panel()])
        self.summary("P3", "KIN3", [panel()])
        self.run_aggregate(panel_min=3)
        reads = self.readouts()
        self.assertEqual(set(reads), {"TR-FRET / HTRF"})                         # the 3-protein panel is left out
        self.assertIn("selectivity panels", reads["TR-FRET / HTRF"]["notes"])
        top = [r for r in self.usage() if r["entry_id"] == "K1" and r["scope_type"] == "gene" and r["scope_id"] == "KIN1"
               and r["effect"] == "inhibition" and not r["readout"]][0]
        self.assertEqual(top["n_campaigns_using"], "2")                          # function occurrence still counts it
        self.run_aggregate(panel_min=0)
        self.assertEqual(set(self.readouts()), {"TR-FRET / HTRF", "Radiometric"})

    def test_old_summary_is_upgraded_once_from_label_file(self):
        old = camp("doc:D1", f=["Kinase activity", "Ligand binding affinity"],
                   t=["TR-FRET / HTRF", "KINOMEscan competition binding"])
        self.summary("P1", "KIN1", [old])
        path = os.path.join(bd.CACHE, "m2_data", "KIN1_labeled.csv")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        rows = [rec("a1", "doc:D1", "Kinase activity", "TR-FRET / HTRF"),
                rec("a2", "doc:D1", "Ligand binding affinity", "KINOMEscan competition binding")]
        with open(path, "w", newline="", encoding="utf-8") as handle:
            w = csv.DictWriter(handle, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        res = self.run_aggregate()
        self.assertEqual(res["upgraded"], 1)
        with open(os.path.join(bd.SUMMARY, "P1.json"), encoding="utf-8") as handle:
            self.assertIn("tpf", json.load(handle)["campaigns"][0])
        self.assertEqual(set(self.readouts()), {"TR-FRET / HTRF"})
        self.assertEqual(self.run_aggregate()["upgraded"], 0)                    # nothing left to upgrade


if __name__ == "__main__":
    unittest.main()
