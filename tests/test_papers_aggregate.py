"""Tests for readout identity in tools/papers_aggregate.py (no network, never touches data/).

A readouts cell is curated free text; a proposal is a controlled label. Comparing them literally
proposed 7 readouts in 27 that the entry already listed (measured on review/11_papers_batches_0_4),
because word order, punctuation, a trailing example and house abbreviations all defeat an exact
match. These tests pin both directions: the duplicates must be recognised, and readouts that only
look similar must still be proposed.
"""
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import papers_aggregate as pa  # noqa: E402


def listed(*readouts):
    return [(r, pa.readout_key(r)) for r in readouts]


class TestDuplicatesAreRecognised(unittest.TestCase):
    """Each case below was a real false proposal in chunk 11."""

    def test_word_order(self):
        self.assertIsNotNone(pa.already_listed("TR-FRET / HTRF", listed("HTRF / TR-FRET cAMP")))
        self.assertIsNotNone(pa.already_listed("HTRF / TR-FRET", listed("TR-FRET / HTRF")))

    def test_trailing_example_is_not_part_of_the_identity(self):
        self.assertIsNotNone(pa.already_listed(
            "Fluorescence polarization",
            listed("fluorescence polarization with a labeled peptide or DNA duplex")))

    def test_house_abbreviation(self):
        # "FP" in the cell, spelled out in the proposal
        self.assertIsNotNone(pa.already_listed(
            "Fluorescence polarization",
            listed("TR-FRET or FP tracer displacement (for example labeled thalidomide for CRBN)")))

    def test_brand_name_means_the_same_measurement(self):
        # Kinase-Glo reads the ATP left unconsumed
        self.assertIsNotNone(pa.already_listed(
            "Luminescent ATP-depletion detection", listed("Kinase-Glo luminescence")))

    def test_implementation_of_a_listed_readout(self):
        self.assertIsNotNone(pa.already_listed("BRET / NanoBRET", listed("BRET")))

    def test_specific_single_token(self):
        self.assertIsNotNone(pa.already_listed(
            "Radioligand binding", listed("Radioligand filtration or scintillation proximity")))


class TestGenuineGapsSurvive(unittest.TestCase):
    """The expensive failure: suppressing a readout the entry really does not have."""

    def test_weak_single_token_does_not_match(self):
        # `luciferase` alone must not match `split-luciferase complementation`: a reporter-gene
        # readout and a split-reporter interaction assay are different assays.
        self.assertIsNone(pa.already_listed(
            "Luciferase (firefly or NanoLuc)", listed("split-luciferase complementation")))

    def test_shared_vendor_suffix_does_not_match(self):
        self.assertIsNone(pa.already_listed("Luminescent Glo assay", listed("ADP-Glo luminescence")))

    def test_unrelated_readouts(self):
        cell = listed("ELISA / immunosorbent", "AlphaScreen / AlphaLISA", "TR-FRET / HTRF")
        for proposed in ("SPR / BLI", "pull-down", "Flow cytometry", "ITC"):
            self.assertIsNone(pa.already_listed(proposed, cell), proposed)

    def test_same_family_different_platform(self):
        # ADP-Glo is not answered by a fluorogenic-substrate cell
        self.assertIsNone(pa.already_listed(
            "ADP-Glo luminescence",
            listed("fluorogenic or luminogenic substrates (Amplex Red, MAO-Glo)")))

    def test_empty_key_never_matches(self):
        # a label made only of stopwords must not swallow every readout
        self.assertEqual(pa.readout_key("binding assay"), frozenset())
        self.assertIsNone(pa.already_listed("binding assay", listed("SPR / BLI")))


class TestReadoutKey(unittest.TestCase):

    def test_order_and_punctuation_do_not_matter(self):
        self.assertEqual(pa.readout_key("TR-FRET / HTRF"), pa.readout_key("HTRF/TR-FRET"))

    def test_parentheses_are_dropped(self):
        self.assertEqual(pa.readout_key("Thermal shift (DSF)"), pa.readout_key("thermal shift"))


if __name__ == "__main__":
    unittest.main()
