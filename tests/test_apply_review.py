"""Tests for tools/apply_review.py on small synthetic files (no network, never touches data/)."""
import contextlib
import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools"))
import apply_review as ar  # noqa: E402

BOM, EOL = "﻿", "\r\n"
ASSAYS = ["entry_id,assay_name,plate_formats,rationale", "X1,Assay one,,old text", 'X2,"Assay, two",96,keep']
KITS = ["kit_id,vendor,product,applies_to_entries,url", "C1,V,Kit one,X1,https://a"]
ENTRY_P = ["entry_id,field,current_value,proposed_value,evidence,confidence,decision,your_note",
           "X1,plate_formats,,96; 384,e,high,accept,",
           "X1,rationale,old text,new text,e,low,reject,",
           "X2,rationale,keep,proposed,e,low,edit,my wording",
           "X2,plate_formats,384,1536,e,low,accept,",          # file holds 96, not 384 -> conflict
           "X1,assay_name,Assay one,Renamed,e,low,,"]           # no decision yet
KIT_P = ["kit_id,vendor,product,applies_to_entries,url,what_claude_read,confidence,decision,your_note",
         "C2,V2,Kit two,X1; X2,https://b,page,high,accept,",
         "C3,V3,Kit three,X9,https://c,page,high,accept,",    # unknown entry -> error
         "C4,V4,Kit four,X1,https://d,page,high,reject,"]


def write(path, lines, bom=True):
    with open(path, "wb") as f:
        f.write(((BOM if bom else "") + EOL.join(lines) + EOL).encode("utf-8"))


class ApplyReviewTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.data = os.path.join(self.tmp.name, "data")
        self.chunk = os.path.join(self.tmp.name, "chunk")
        os.makedirs(self.data); os.makedirs(self.chunk)
        write(os.path.join(self.data, "assays.csv"), ASSAYS)
        write(os.path.join(self.data, "kits.csv"), KITS)
        write(os.path.join(self.chunk, "entry_proposals.csv"), ENTRY_P, bom=False)
        write(os.path.join(self.chunk, "kit_proposals.csv"), KIT_P, bom=False)

    def tearDown(self):
        self.tmp.cleanup()

    def run_tool(self, *extra):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            ar.main([self.chunk, "--data", self.data] + list(extra))
        return out.getvalue()

    def raw(self, name):
        with open(os.path.join(self.data, name), "rb") as f:
            return f.read()

    def test_dry_run_writes_nothing(self):
        before = self.raw("assays.csv"), self.raw("kits.csv")
        out = self.run_tool()
        self.assertEqual(before, (self.raw("assays.csv"), self.raw("kits.csv")))
        self.assertIn("CONFLICT", out)
        self.assertIn("ERROR", out)
        self.assertFalse(os.path.exists(os.path.join(self.data, "backup")))

    def test_apply_changes_only_accepted_rows_and_keeps_layout(self):
        self.run_tool("--apply")
        _, rows, lay = ar.read_table(os.path.join(self.data, "assays.csv"))
        r = {a["entry_id"]: a for a in rows}
        self.assertEqual(r["X1"]["plate_formats"], "96; 384")      # accepted
        self.assertEqual(r["X1"]["rationale"], "old text")         # rejected
        self.assertEqual(r["X2"]["rationale"], "my wording")       # edit uses your_note
        self.assertEqual(r["X2"]["plate_formats"], "96")           # conflict left alone
        self.assertEqual(r["X1"]["assay_name"], "Assay one")       # undecided left alone
        self.assertEqual(lay, {"bom": True, "eol": "\r\n"})
        self.assertIn(b'"Assay, two"', self.raw("assays.csv"))
        _, kits, _ = ar.read_table(os.path.join(self.data, "kits.csv"))
        self.assertEqual([k["kit_id"] for k in kits], ["C1", "C2"])  # C3 error, C4 rejected
        self.assertEqual(set(kits[1].keys()), {"kit_id", "vendor", "product", "applies_to_entries", "url"})
        self.assertEqual(len(os.listdir(os.path.join(self.data, "backup"))), 2)

    def test_second_run_is_a_no_op(self):
        self.run_tool("--apply")
        after = self.raw("assays.csv"), self.raw("kits.csv")
        out = self.run_tool("--apply")
        self.assertEqual(after, (self.raw("assays.csv"), self.raw("kits.csv")))
        self.assertIn("Nothing to write", out)

    def test_kit_updates_change_existing_kits(self):
        write(os.path.join(self.chunk, "kit_updates.csv"), [
            "kit_id,field,current_value,proposed_value,evidence,confidence,decision,your_note",
            "C1,applies_to_entries,X1,X1; X2,e,high,accept,",
            "C1,url,https://OLD,https://new,e,high,accept,",       # stale current_value -> conflict
            "C1,vendor,V,V9,e,high,accept,",
            "C9,url,,https://x,e,high,accept,"], bom=False)       # unknown kit -> error
        out = self.run_tool("--apply")
        self.assertIn("CONFLICT", out)
        _, kits, _ = ar.read_table(os.path.join(self.data, "kits.csv"))
        k = {r["kit_id"]: r for r in kits}
        self.assertEqual(k["C1"]["applies_to_entries"], "X1; X2")
        self.assertEqual(k["C1"]["vendor"], "V9")
        self.assertEqual(k["C1"]["url"], "https://a")
        self.assertEqual(sorted(k), ["C1", "C2"])

    def test_kit_update_to_unknown_entry_is_refused(self):
        write(os.path.join(self.chunk, "kit_updates.csv"), [
            "kit_id,field,current_value,proposed_value,evidence,confidence,decision,your_note",
            "C1,applies_to_entries,X1,X1; NOPE,e,high,accept,"], bom=False)
        self.run_tool("--apply")
        _, kits, _ = ar.read_table(os.path.join(self.data, "kits.csv"))
        self.assertEqual(kits[0]["applies_to_entries"], "X1")

    def test_new_classes_and_entries_and_class_updates(self):
        write(os.path.join(self.data, "function_classes.csv"),
              ["class_id,label,keywords", "kin,Kinase,kinase", "generic,Any,"])
        write(os.path.join(self.chunk, "new_classes.csv"), [
            "class_id,label,keywords,evidence,decision,your_note",
            "parp,PARP,ADP-ribosyl,e,accept,"], bom=False)
        write(os.path.join(self.chunk, "class_updates.csv"), [
            "class_id,field,current_value,proposed_value,evidence,confidence,decision,your_note",
            "kin,keywords,kinase,Serine/threonine-protein kinase,e,high,accept,"], bom=False)
        write(os.path.join(self.chunk, "new_entries.csv"), [
            "entry_id,assay_name,plate_formats,rationale,classes,decision,your_note",
            "PR1,PARP activity,,new,parp,accept,",           # uses the class added in the same chunk
            "PR9,Bad entry,,x,nosuchclass,accept,"], bom=False)
        # the fake assays.csv has no classes column, so give it one for this test
        write(os.path.join(self.data, "assays.csv"),
              ["entry_id,assay_name,plate_formats,rationale,classes", "X1,Assay one,,old text,kin", 'X2,"Assay, two",96,keep,kin'])
        out = self.run_tool("--apply")
        self.assertIn("uses unknown classes: nosuchclass", out)
        _, classes, _ = ar.read_table(os.path.join(self.data, "function_classes.csv"))
        self.assertEqual([c["class_id"] for c in classes], ["kin", "generic", "parp"])
        self.assertEqual(classes[0]["keywords"], "Serine/threonine-protein kinase")
        _, assays, _ = ar.read_table(os.path.join(self.data, "assays.csv"))
        self.assertIn("PR1", [a["entry_id"] for a in assays])
        self.assertNotIn("PR9", [a["entry_id"] for a in assays])

    def test_allowed_new_column_is_added_but_unknown_column_is_refused(self):
        write(os.path.join(self.chunk, "entry_proposals.csv"), [
            "entry_id,field,current_value,proposed_value,evidence,confidence,decision,your_note",
            "X1,qc_metric,,hts,AGM,high,accept,",
            "X2,dmso_pct,,1-2,AGM,high,accept,",
            "X1,made_up_column,,x,e,high,accept,"], bom=False)       # not in NEW_COLUMNS -> error
        out = self.run_tool("--apply")
        self.assertIn("ERROR", out)
        h, rows, lay = ar.read_table(os.path.join(self.data, "assays.csv"))
        self.assertEqual(h[-2:], ["qc_metric", "dmso_pct"])
        self.assertNotIn("interference", h)                         # only columns a proposal uses are added
        self.assertNotIn("made_up_column", h)
        r = {a["entry_id"]: a for a in rows}
        self.assertEqual((r["X1"]["qc_metric"], r["X1"]["dmso_pct"]), ("hts", ""))
        self.assertEqual((r["X2"]["qc_metric"], r["X2"]["dmso_pct"]), ("", "1-2"))
        self.assertEqual(lay, {"bom": True, "eol": "\r\n"})
        self.assertIn("Nothing to write", self.run_tool("--apply"))  # second run: already applied

    def test_unchanged_file_round_trips_byte_for_byte(self):
        path = os.path.join(self.data, "assays.csv")
        before = self.raw("assays.csv")
        h, rows, lay = ar.read_table(path)
        ar.write_table(path, h, rows, lay)
        self.assertEqual(before, self.raw("assays.csv"))


if __name__ == "__main__":
    unittest.main()
