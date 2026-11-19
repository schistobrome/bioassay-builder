"""Tests for the term_proposals.csv -> scripts/term_table.csv route in tools/apply_review.py.

Vocabulary was the one finding of the paper crawl with no way into the repository: a term could be
proposed and decided but never written anywhere. This route closes that, and these tests pin the
parts that make it safe to use on the labeller's own file -- a term is written only on an explicit
accept, a term that cannot be filed is refused by name rather than guessed at, and running it twice
does nothing the second time.

Every fixture is self-contained: --data points at a temporary copy, so the real
scripts/term_table.csv is never touched.
"""
import contextlib
import io
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import apply_review as ar  # noqa: E402

BOM, EOL = "﻿", "\r\n"

TERM_TABLE = ["axis,term,synonyms,examples,suppresses,applies_to,notes",
              "technology,ADP-Glo luminescence,adpglo,ADP-Glo readout,,any,"]

COLS = ("axis,term,synonyms,examples,suppresses,applies_to,notes,"
        "term_id,kind,n_papers,evidence,confidence,decision,your_note")
TERM_P = [COLS,
          # accepted: a straightforward new technology
          "technology,Microscale thermophoresis (MST),mst|thermophoresis,MST binding,,any,note,"
          "T003,technology,10,q,high,accept,",
          # rejected
          "technology,Alkaline comet assay,comet,comet tails,,any,note,T076,technology,2,q,medium,reject,",
          # no decision
          "function,Protein lysine lactylation,lactylation,lactyl,,any,note,F011,function,2,q,medium,,",
          # a system term: no axis, must be refused by name, not filed somewhere plausible
          ",Ba/F3 cytokine-independence model,baf3,baf3 cells,,any,note,S001,system,3,q,high,accept,",
          # an axis m2_label.py would never match on
          "preparation,Crude membrane vesicles,vesicles,vesicle prep,,any,note,S009,system,2,q,high,accept,",
          # already present in the table -> idempotent, not a duplicate row
          "technology,ADP-Glo luminescence,adpglo,ADP-Glo readout,,any,,T999,technology,5,q,high,accept,"]


def write(path, lines, bom=True):
    with open(path, "wb") as f:
        f.write(((BOM if bom else "") + EOL.join(lines) + EOL).encode("utf-8"))


def read(path):
    with open(path, encoding="utf-8-sig") as f:
        return f.read()


class TermRouteTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.data = os.path.join(self.tmp.name, "data")
        self.chunk = os.path.join(self.tmp.name, "chunk")
        os.makedirs(self.data)
        os.makedirs(self.chunk)
        # a local term_table.csv, so OTHER_DIRS never reaches the real scripts/ copy
        write(os.path.join(self.data, "term_table.csv"), TERM_TABLE)
        write(os.path.join(self.chunk, "term_proposals.csv"), TERM_P, bom=False)
        self.table = os.path.join(self.data, "term_table.csv")

    def tearDown(self):
        self.tmp.cleanup()

    def run_tool(self, *extra):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = ar.main([self.chunk, "--data", self.data] + list(extra))
        return rc, buf.getvalue()

    def test_dry_run_writes_nothing(self):
        before = read(self.table)
        _, out = self.run_tool()
        self.assertIn("term_proposals.csv -> ", out)
        self.assertEqual(read(self.table), before)

    def test_accepted_term_is_written(self):
        self.run_tool("--apply")
        body = read(self.table)
        self.assertIn("Microscale thermophoresis (MST)", body)
        self.assertIn("mst|thermophoresis", body)

    def test_rejected_and_undecided_terms_are_not_written(self):
        self.run_tool("--apply")
        body = read(self.table)
        self.assertNotIn("Alkaline comet", body)
        self.assertNotIn("lactylation", body)

    def test_term_with_no_axis_is_refused_by_name(self):
        _, out = self.run_tool("--apply")
        self.assertIn("Ba/F3", out)
        self.assertIn("no axis", out)
        self.assertNotIn("Ba/F3", read(self.table))

    def test_unknown_axis_is_refused(self):
        _, out = self.run_tool("--apply")
        self.assertIn("axis must be one of", out)
        self.assertNotIn("Crude membrane vesicles", read(self.table))

    def test_existing_term_is_not_duplicated(self):
        self.run_tool("--apply")
        self.assertEqual(read(self.table).count("ADP-Glo luminescence"), 1)

    def test_running_twice_changes_nothing(self):
        self.run_tool("--apply")
        after_first = read(self.table)
        _, out = self.run_tool("--apply")
        self.assertEqual(read(self.table), after_first)
        self.assertIn("Nothing to write.", out)

    def test_the_relabel_warning_is_printed(self):
        _, out = self.run_tool("--apply")
        self.assertIn("relabel", out.lower())

    def test_a_chunk_from_an_older_version_is_reported_not_crashed(self):
        """term_proposals.csv gained the term_table.csv columns; chunks rendered before that do
        not have them. The tool must say so and leave the file alone, not raise KeyError."""
        old = ["term_id,proposed_term,kind,definition,variants,n_papers,evidence,confidence,"
               "decision,your_note",
               "T003,Microscale thermophoresis (MST),technology,d,mst,10,q,high,accept,"]
        write(os.path.join(self.chunk, "term_proposals.csv"), old, bom=False)
        before = read(self.table)
        rc, out = self.run_tool("--apply")
        self.assertIn("older version", out)
        self.assertIn("re-render", out)
        self.assertEqual(read(self.table), before)


if __name__ == "__main__":
    unittest.main()
