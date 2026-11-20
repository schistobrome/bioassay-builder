"""Tests for the readout-label check in tools/build_db.py (no network, never touches data/).

The page matches an assay's curated readouts to usage readout labels with a regex per label. That
list used to be hard-coded in docs/index.html and kept in step with scripts/term_table.csv by hand:
32 terms, 32 entries, nothing enforcing it. Accepting a technology term and re-running the relabel
would put a label into the data that the page could not match, silently. These tests pin the
check that now stands between the two.
"""
import csv
import os
import re
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import build_db  # noqa: E402

check = build_db.check_readout_labels


def rows(*pairs):
    return [{"label": l, "pattern": p, "notes": ""} for l, p in pairs]


def load(path):
    with open(path, newline="", encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))


class TestTheDataMustBeMatchable(unittest.TestCase):

    def test_clean_when_all_three_agree(self):
        e, w, out = check(rows(("ITC", r"\bitc\b")), {"ITC": 5}, {"ITC"})
        self.assertEqual((e, w), ([], []))
        self.assertEqual(out, [["ITC", r"\bitc\b"]])

    def test_a_label_in_the_data_with_no_row_stops_the_build(self):
        e, w, _ = check(rows(("ITC", r"\bitc\b")), {"ITC": 5, "Microscale thermophoresis": 12},
                        {"ITC", "Microscale thermophoresis"})
        self.assertEqual(len(e), 1)
        self.assertIn("Microscale thermophoresis", e[0])
        self.assertIn("12 rows", e[0])
        self.assertIn("Add: Microscale thermophoresis,", e[0])
        self.assertEqual(w, [], "the same gap must not be reported twice")

    def test_an_accepted_term_not_yet_in_the_data_only_warns(self):
        """Between accepting a term and re-running the relabel the term is inert: say what to add,
        but do not block every other build in the meantime."""
        e, w, _ = check(rows(("ITC", r"\bitc\b")), {"ITC": 5}, {"ITC", "GloSensor cAMP biosensor"})
        self.assertEqual(e, [])
        self.assertEqual(len(w), 1)
        self.assertIn("GloSensor cAMP biosensor", w[0])
        self.assertIn("before the m2 relabel", w[0])

    def test_a_row_the_vocabulary_no_longer_has_warns(self):
        e, w, _ = check(rows(("ITC", r"\bitc\b"), ("Old name", "old")), {"ITC": 5}, {"ITC"})
        self.assertEqual(e, [])
        self.assertEqual(len(w), 1)
        self.assertIn("Old name", w[0])

    def test_without_term_table_only_the_data_is_checked(self):
        e, w, _ = check(rows(("ITC", r"\bitc\b")), {"ITC": 5}, None)
        self.assertEqual((e, w), ([], []))


class TestTheTableItself(unittest.TestCase):

    def test_duplicate_label(self):
        e, _, out = check(rows(("ITC", "itc"), ("ITC", "calorimetr")), {}, None)
        self.assertTrue(any("duplicate" in x for x in e))
        self.assertEqual(len(out), 1)

    def test_missing_label_or_pattern(self):
        e, _, out = check(rows(("", "x"), ("ITC", "")), {}, None)
        self.assertEqual(len(e), 2)
        self.assertEqual(out, [])

    def test_pattern_that_does_not_compile(self):
        e, _, out = check(rows(("ITC", "(itc")), {}, None)
        self.assertIn("does not compile", e[0])
        self.assertEqual(out, [])

    def test_python_only_syntax_is_refused(self):
        """These compile in Python and fail in the browser, which is where the pattern runs."""
        for pat in ["(?i)itc", "(?P<x>itc)", r"itc\Z", "(?>itc)", "it+c*+"]:
            e, _, out = check(rows(("ITC", pat)), {}, None)
            self.assertTrue(e and "browser" in e[0], pat)
            self.assertEqual(out, [], pat)

    def test_js_lookbehind_is_allowed(self):
        e, _, _ = check(rows(("Thermal shift (DSF)", r"(?<!cellular )thermal[ -]shift")), {}, None)
        self.assertEqual(e, [])

    def test_order_is_kept(self):
        _, _, out = check(rows(("B", "b"), ("A", "a"), ("C", "c")), {}, None)
        self.assertEqual([l for l, _ in out], ["B", "A", "C"])


class TestTheShippedTable(unittest.TestCase):
    """The real data/readout_labels.csv against the real scripts/term_table.csv."""

    @classmethod
    def setUpClass(cls):
        cls.rows = load(os.path.join(ROOT, "data", "readout_labels.csv"))
        cls.terms = {r["term"] for r in load(os.path.join(ROOT, "scripts", "term_table.csv"))
                     if (r.get("axis") or "").strip().lower() == "technology"}
        cls.rx = {r["label"]: re.compile(r["pattern"], re.I) for r in cls.rows}

    def test_one_row_per_technology_term(self):
        e, w, out = check(self.rows, {}, self.terms)
        self.assertEqual(e, [])
        self.assertEqual(w, [], "\n".join(w))
        self.assertEqual({l for l, _ in out}, self.terms)

    def test_moved_patterns_still_discriminate(self):
        """Lookarounds that only work if the move out of the JS literal kept them intact."""
        fret = self.rx["Fluorogenic / FRET substrate"]
        self.assertTrue(fret.search("FRET peptide substrate"))
        self.assertFalse(fret.search("TR-FRET"))
        self.assertFalse(fret.search("FRET / HTRF"))
        dsf = self.rx["Thermal shift (DSF)"]
        self.assertTrue(dsf.search("thermal shift assay"))
        self.assertFalse(dsf.search("cellular thermal shift"))


if __name__ == "__main__":
    unittest.main()
