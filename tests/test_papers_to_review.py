"""Tests for tools/papers_to_review.py on small synthetic files (no network, never touches data/).

The behaviour under test is the one that protects the owner's work: a chunk is regenerated every
time a batch lands, and a regeneration must not renumber a drafted entry or forget a decision
already recorded against it -- not even when the candidate behind it temporarily disappears.
"""
import contextlib
import csv
import io
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import papers_to_review as ptr  # noqa: E402

ASSAYS = ["entry_id,assay_name,classes,assay_category,scope,role,screening_suitability,"
          "biological_relevance,measures,readouts,evidence_labels,evidence_function,conditions,"
          "reagents,rationale,status,development_priority,rating_notes,plate_formats,sources_doi,"
          "learn_more_url,reviewer,reviewed_on",
          "K1,Kinase activity,protein_kinase,biochemical_hts,cell-free,primary,high,high,x,"
          "ADP-Glo luminescence,Kinase activity,Kinase activity,,,,,,,,,,,"]

# Two candidates; the second one is the one we will make disappear and bring back.
ENTRY_CANDIDATES = [
    "class_id,scope,function_label,n_papers,n_proteins,readouts_seen,proteins,pmcids,quote",
    "protein_kinase,cell-based,Target engagement (cellular),13,6,CETSA,AKT1; EGFR,PMC1; PMC2,q",
    "protein_kinase,cell-based,Cellular signaling,6,5,Western blot,AKT1; JAK2,PMC3,q",
]
# Header plus the SECOND candidate only, so the first one (which becomes LIT1) is the one that goes.
WITHOUT_LIT1 = [ENTRY_CANDIDATES[0], ENTRY_CANDIDATES[2]]


def write(path, lines):
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def rows(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


class PapersToReviewTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.data = os.path.join(self.tmp.name, "data")
        self.papers = os.path.join(self.tmp.name, "papers")
        self.chunk = os.path.join(self.tmp.name, "chunk")
        for d in (self.data, self.papers):
            os.makedirs(d)
        write(os.path.join(self.data, "assays.csv"), ASSAYS)
        self.candidates(ENTRY_CANDIDATES)
        self._saved = ptr.DATA, ptr.PAPERS
        ptr.DATA, ptr.PAPERS = self.data, self.papers

    def tearDown(self):
        ptr.DATA, ptr.PAPERS = self._saved
        self.tmp.cleanup()

    def candidates(self, entry_lines):
        write(os.path.join(self.papers, "new_entry_candidates.csv"), entry_lines)
        write(os.path.join(self.papers, "readout_candidates.csv"),
              ["class_id,scope,readout,n_papers,n_proteins,entries_seen,proteins,pmcids,quote",
               "protein_kinase,cell-free,SPR / BLI,15,7,K1,AKT1,PMC9,q"])
        write(os.path.join(self.papers, "term_candidates.csv"),
              ["term_id,proposed_term,kind,nearest_controlled_value,definition,variants,n_papers,"
               "n_proteins,proteins,pmcids,quote,catalog_status",
               "T1,DARTS,technology,CETSA,d,v,2,2,AKT1,PMC1,q,supported"])

    def run_tool(self):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            ptr.main(["--chunk", self.chunk])
        return out.getvalue()

    def entry(self, eid):
        for r in rows(os.path.join(self.chunk, "new_entries.csv")):
            if r["entry_id"] == eid:
                return r
        return None

    def decide(self, eid, **fields):
        path = os.path.join(self.chunk, "new_entries.csv")
        data = rows(path)
        header = list(data[0].keys())
        for r in data:
            if r["entry_id"] == eid:
                r.update(fields)
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=header)
            w.writeheader()
            w.writerows(data)

    def test_ids_are_keyed_to_the_entry_not_its_position(self):
        self.run_tool()
        first = {r["entry_id"]: r["evidence_labels"]
                 for r in rows(os.path.join(self.chunk, "new_entries.csv"))}
        # Reverse the candidate order: the same entries must keep the same ids.
        self.candidates([ENTRY_CANDIDATES[0], ENTRY_CANDIDATES[2], ENTRY_CANDIDATES[1]])
        self.run_tool()
        again = {r["entry_id"]: r["evidence_labels"]
                 for r in rows(os.path.join(self.chunk, "new_entries.csv"))}
        self.assertEqual(first, again)

    def test_a_decision_survives_regeneration(self):
        self.run_tool()
        self.decide("LIT1", decision="accept", screening_suitability="moderate",
                    biological_relevance="high", assay_name="My own wording")
        out = self.run_tool()
        r = self.entry("LIT1")
        self.assertEqual(r["decision"], "accept")
        self.assertEqual(r["screening_suitability"], "moderate")
        self.assertEqual(r["assay_name"], "My own wording")
        self.assertIn("kept", out)

    def test_a_decision_survives_the_candidate_disappearing_and_coming_back(self):
        """The hazard this exists for: a narrower --batches or a stricter threshold drops a
        candidate, and the id it held must not be handed to a different function later."""
        self.run_tool()
        key_before = self.entry("LIT1")["evidence_labels"]
        self.decide("LIT1", decision="accept", screening_suitability="low",
                    biological_relevance="high")
        self.candidates(WITHOUT_LIT1)           # LIT1's candidate is gone
        self.run_tool()
        self.assertIsNone(self.entry("LIT1"))
        self.candidates(ENTRY_CANDIDATES)       # and back
        self.run_tool()
        r = self.entry("LIT1")
        self.assertEqual(r["evidence_labels"], key_before)
        self.assertEqual(r["decision"], "accept")
        self.assertEqual(r["screening_suitability"], "low")

    def test_evidence_refreshes_while_the_decision_is_kept(self):
        """New papers must still update the counts: only the owner's own fields are frozen."""
        self.run_tool()
        self.decide("LIT1", decision="accept")
        grown = list(ENTRY_CANDIDATES)
        grown[1] = ("protein_kinase,cell-based,Target engagement (cellular),99,9,CETSA,"
                    "AKT1; EGFR,PMC1; PMC2,q")
        self.candidates(grown)
        self.run_tool()
        r = self.entry("LIT1")
        self.assertEqual(r["decision"], "accept")
        self.assertIn("99", r["rationale"])

    def test_ratings_are_blank_until_a_person_fills_them(self):
        self.run_tool()
        r = self.entry("LIT1")
        self.assertEqual(r["screening_suitability"], "")
        self.assertEqual(r["biological_relevance"], "")
        self.assertEqual(r["assay_category"], "biophysical_binding")   # inferred, not blank
        self.assertEqual(r["role"], "confirmation")


if __name__ == "__main__":
    unittest.main()
