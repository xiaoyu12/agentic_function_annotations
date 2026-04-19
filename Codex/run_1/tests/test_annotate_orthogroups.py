import unittest
from pathlib import Path

import annotate_orthogroups as ao


ROOT = Path(__file__).resolve().parents[1]


class AnnotationParserTests(unittest.TestCase):
    def test_load_signalp_predictions_reads_positive_call(self):
        predictions = ao.load_signalp_predictions(ROOT / "signalp_results")

        rec = predictions["OG0000049__000328"]
        self.assertEqual(rec["prediction"], "Signal Peptide (Sec/SPI)")
        self.assertAlmostEqual(rec["sp_probability"], 0.9848, places=3)
        self.assertEqual(rec["cleavage_site_end"], 26)

    def test_load_deeptmhmm_predictions_reads_tm_count(self):
        predictions = ao.load_deeptmhmm_predictions(ROOT / "deeptmhmm_results")

        rec = predictions["OG0000049__000112"]
        self.assertEqual(rec["class"], "TM")
        self.assertEqual(rec["tm_count"], 6)

    def test_parse_pfam_tbl_extracts_domain_names(self):
        hits = ao.parse_pfam_tbl(ROOT / "hmm_out" / "OG0000049.pfam.tbl")

        domain_names = {hit["domain_name"] for hit in hits}
        self.assertIn("FKBP15", domain_names)
        self.assertIn("CCDC39", domain_names)

    def test_filter_pfam_hits_drops_weak_reporting_only_hits(self):
        hits = ao.parse_pfam_tbl(ROOT / "hmm_out" / "OG0009246.pfam.tbl")

        filtered = ao.filter_pfam_hits(hits)
        self.assertEqual(filtered, [])

    def test_parse_blast_tsv_extracts_subject_description(self):
        hits = ao.parse_blast_tsv(ROOT / "blast_out" / "OG0009301.blast.tsv")

        self.assertIn("DNA-directed RNA polymerase II subunit RPB1", hits[0]["subject_description"])

    def test_filter_reliable_blast_hits_drops_weak_single_hit(self):
        hits = ao.parse_blast_tsv(ROOT / "blast_out" / "OG0009301.blast.tsv")

        filtered = ao.filter_reliable_blast_hits(hits)
        self.assertEqual(filtered, [])


if __name__ == "__main__":
    unittest.main()
