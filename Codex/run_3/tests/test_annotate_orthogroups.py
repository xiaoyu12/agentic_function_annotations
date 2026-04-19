import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "annotate_orthogroups.py"


def load_module(testcase: unittest.TestCase):
    testcase.assertTrue(
        MODULE_PATH.exists(),
        f"Expected annotation module at {MODULE_PATH}",
    )
    spec = importlib.util.spec_from_file_location("annotate_orthogroups", MODULE_PATH)
    testcase.assertIsNotNone(spec)
    testcase.assertIsNotNone(spec.loader)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AnnotateOrthogroupsTests(unittest.TestCase):
    def test_parse_blast_file_summarizes_clean_descriptions(self):
        module = load_module(self)

        summary = module.parse_blast_file(ROOT / "blast_out" / "OG0001332.blast.tsv")

        self.assertEqual(summary["hit_count"], 191)
        self.assertEqual(summary["top_descriptions"][0]["description"], "Protein hedgehog")
        self.assertEqual(summary["top_descriptions"][0]["count"], 81)

    def test_parse_pfam_tbl_aggregates_related_domain_hits(self):
        module = load_module(self)

        summary = module.parse_pfam_tbl(ROOT / "hmm_out" / "OG0020700.pfam.tbl")

        self.assertEqual(summary["domain_hit_count"], 7)
        self.assertEqual(summary["top_domains"][0]["description"], "PDZ domain")
        self.assertEqual(summary["top_domains"][0]["count"], 7)

    def test_topology_and_signalp_results_are_grouped_by_orthogroup(self):
        module = load_module(self)

        signalp = module.load_signalp_results(ROOT / "signalp_results")
        deeptmhmm = module.load_deeptmhmm_results(ROOT / "deeptmhmm_results")

        self.assertEqual(signalp["OG0009301"]["positive_count"], 21)
        self.assertEqual(signalp["OG0009301"]["total_count"], 29)
        self.assertEqual(deeptmhmm["OG0009301"]["type_counts"]["SP+TM"], 11)
        self.assertEqual(deeptmhmm["OG0009301"]["type_counts"]["SP"], 11)
        self.assertEqual(deeptmhmm["OG0009301"]["type_counts"]["TM"], 4)
        self.assertEqual(deeptmhmm["OG0009301"]["type_counts"]["GLOB"], 3)
        self.assertEqual(deeptmhmm["OG0009301"]["tm_sequence_count"], 15)
        self.assertEqual(deeptmhmm["OG0009301"]["signal_peptide_count"], 22)

    def test_build_evidence_rows_covers_all_orthogroups(self):
        module = load_module(self)

        rows = module.build_evidence_rows(ROOT)
        row = next(item for item in rows if item["orthogroup"] == "OG0014155")

        self.assertEqual(len(rows), 73)
        self.assertEqual(row["sequence_count"], 22)
        self.assertEqual(
            row["top_blast_description"],
            "Pyrimidine-specific ribonucleoside hydrolase RihA",
        )
        self.assertEqual(
            row["top_pfam_description"],
            "Inosine-uridine preferring nucleoside hydrolase",
        )
        self.assertEqual(row["signalp_positive_count"], 11)

    def test_load_manual_annotations_covers_all_orthogroups(self):
        module = load_module(self)

        annotations = module.load_manual_annotations(
            ROOT / "annotations" / "orthogroup_manual_annotations.tsv"
        )

        self.assertEqual(len(annotations), 73)
        self.assertEqual(
            annotations["OG0009246"]["predicted_function"],
            "Secreted adhesive ECM/collagen-like protein",
        )
        self.assertEqual(annotations["OG0009246"]["calcification_relevance"], "high")

    def test_build_annotated_rows_merges_curated_calls(self):
        module = load_module(self)

        rows = module.build_annotated_rows(
            ROOT, ROOT / "annotations" / "orthogroup_manual_annotations.tsv"
        )
        row = next(item for item in rows if item["orthogroup"] == "OG0023594")

        self.assertEqual(
            row["predicted_function"],
            "GT8-family glycosyltransferase with ankyrin-repeat scaffold features",
        )
        self.assertEqual(row["calcification_relevance"], "high")
        self.assertIn("matrix polysaccharides", row["calcification_rationale"])

    def test_write_outputs_creates_summary_tsv_and_markdown_report(self):
        module = load_module(self)

        rows = module.build_annotated_rows(
            ROOT, ROOT / "annotations" / "orthogroup_manual_annotations.tsv"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            module.write_outputs(rows, output_dir)

            summary_path = output_dir / "orthogroup_annotation_summary.tsv"
            report_path = output_dir / "orthogroup_annotation_report.md"

            self.assertTrue(summary_path.exists())
            self.assertTrue(report_path.exists())

            report_text = report_path.read_text()
            self.assertIn("Calcification-Relevant Orthogroup Annotations", report_text)
            self.assertIn("High-priority calcification candidates", report_text)
            self.assertIn("OG0009246", report_text)


if __name__ == "__main__":
    unittest.main()
