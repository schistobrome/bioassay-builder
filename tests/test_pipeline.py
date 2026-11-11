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


if __name__ == "__main__":
    unittest.main()
