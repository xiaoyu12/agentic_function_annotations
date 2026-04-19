import tempfile
import textwrap
import unittest
from pathlib import Path

from orthogroup_annotation_report import (
    infer_calcification_relevance,
    infer_function,
    orthogroup_from_surrogate_id,
    normalize_blast_description,
    parse_deeptmhmm_3line,
    parse_pfam_tbl,
)


class OrthogroupAnnotationReportTests(unittest.TestCase):
    def test_orthogroup_from_surrogate_id_uses_prefix_before_double_underscore(self):
        self.assertEqual(orthogroup_from_surrogate_id("OG0013037__000003"), "OG0013037")

    def test_normalize_blast_description_strips_accession_and_species_metadata(self):
        raw = (
            "sp|Q8IDX6|RH2A_PLAF7 Reticulocyte-binding protein homolog 2a "
            "OS=Plasmodium falciparum (isolate 3D7) OX=36329 GN=RH2a"
        )
        self.assertEqual(
            normalize_blast_description(raw),
            "Reticulocyte-binding protein homolog 2a",
        )

    def test_parse_pfam_tbl_reads_domain_hits(self):
        content = textwrap.dedent(
            """\
            # target name accession query name accession E-value score bias E-value score bias exp reg clu ov env dom rep inc description of target
            FKBP15 PF23649.2 jgi|Pleelo874_1|6485806|MIX16057_2_5 - 3.4e-05 24.2 19.6 3.4e-05 24.2 19.6 3.2 3 1 1 4 4 1 1 FK506-binding protein 15-like domain
            """
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sample.pfam.tbl"
            path.write_text(content)
            rows = parse_pfam_tbl(path)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["query"], "jgi|Pleelo874_1|6485806|MIX16057_2_5")
        self.assertEqual(rows[0]["domain"], "FKBP15")
        self.assertEqual(rows[0]["description"], "FK506-binding protein 15-like domain")

    def test_parse_deeptmhmm_3line_reads_label_and_sequence_id(self):
        content = textwrap.dedent(
            """\
            >OG0000049__000001 | GLOB
            MSGIDSPSEQERERRL
            IIIIIIIIIIIIIIII
            >OG0000049__000002 | SP+TM
            MAAAAAAVVVVV
            MMMMSSSSSOOOO
            """
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sample.predicted_topologies.3line"
            path.write_text(content)
            rows = parse_deeptmhmm_3line(path)
        self.assertEqual(
            rows,
            [
                {
                    "sequence_id": "OG0000049__000001",
                    "label": "GLOB",
                    "length": 16,
                    "topology": "IIIIIIIIIIIIIIII",
                },
                {
                    "sequence_id": "OG0000049__000002",
                    "label": "SP+TM",
                    "length": 12,
                    "topology": "MMMMSSSSSOOOO",
                },
            ],
        )

    def test_infer_function_prefers_supported_pfam_family_over_noisy_blast_hit(self):
        row = {
            "top_pfam_description": "FK506-binding protein 15-like domain",
            "top_blast_description": "Reticulocyte-binding protein homolog 2a",
            "top_pfam_count": 7,
            "pfam_annotated_queries": 14,
            "top_blast_count": 6,
            "blast_annotated_queries": 20,
            "signalp_fraction": 0.004,
            "tm_fraction": 0.006,
            "glob_fraction": 0.986,
            "top_pfam_summary": "FK506-binding protein 15-like domain (7)",
            "top_blast_summary": "Reticulocyte-binding protein homolog 2a (6)",
        }
        self.assertEqual(infer_function(row), "FK506-binding protein 15-like protein")

    def test_infer_calcification_relevance_marks_glycosyltransferase_as_plausible_indirect(self):
        row = {
            "function_call": "GT47 exostosin-like membrane glycosyltransferase",
            "top_pfam_description": "Exostosin GT47 domain",
            "top_blast_description": "",
            "top_pfam_summary": "Exostosin GT47 domain (4)",
            "top_blast_summary": "",
            "signalp_fraction": 0.333,
            "tm_fraction": 1.0,
        }
        relevance, rationale = infer_calcification_relevance(row)
        self.assertEqual(relevance, "plausible indirect")
        self.assertIn("trafficking", rationale.lower())

    def test_infer_function_prefers_strong_blast_when_pfam_only_gives_generic_domain(self):
        row = {
            "top_pfam_description": "AAA domain",
            "top_blast_description": "ATP-dependent RecD2 DNA helicase",
            "top_pfam_count": 10,
            "pfam_annotated_queries": 11,
            "top_blast_count": 11,
            "blast_annotated_queries": 11,
            "signalp_fraction": 0.091,
            "tm_fraction": 0.0,
            "glob_fraction": 1.0,
            "top_pfam_summary": "AAA domain (10); Viral superfamily 1 RNA helicase core domain (10)",
            "top_blast_summary": "ATP-dependent RecD2 DNA helicase (11)",
        }
        self.assertEqual(infer_function(row), "Likely ATP-dependent RecD2 DNA helicase")

    def test_infer_function_drops_weak_nonmembrane_blast_name_for_tm_rich_family(self):
        row = {
            "top_pfam_description": "",
            "top_blast_description": "Hexamerin 110",
            "top_pfam_count": 0,
            "pfam_annotated_queries": 0,
            "top_blast_count": 3,
            "blast_annotated_queries": 6,
            "signalp_fraction": 0.056,
            "tm_fraction": 0.833,
            "glob_fraction": 0.0,
            "top_pfam_summary": "",
            "top_blast_summary": "Hexamerin 110 (3); Collagen alpha-2(IV) chain (1)",
        }
        self.assertEqual(infer_function(row), "Membrane protein family with no strong conserved annotation")

    def test_infer_function_uses_sparse_but_concordant_sulfotransferase_evidence(self):
        row = {
            "top_pfam_description": "Sulfotransferase family",
            "top_blast_description": "Heparan-sulfate 6-O-sulfotransferase 1",
            "top_pfam_count": 1,
            "pfam_annotated_queries": 1,
            "top_blast_count": 2,
            "blast_annotated_queries": 2,
            "signalp_fraction": 0.364,
            "tm_fraction": 0.0,
            "glob_fraction": 0.636,
            "top_pfam_summary": "Sulfotransferase family (1); Sulfotransferase domain (1)",
            "top_blast_summary": "Heparan-sulfate 6-O-sulfotransferase 1 (2)",
        }
        self.assertEqual(infer_function(row), "heparan-sulfate sulfotransferase-like enzyme")

    def test_infer_calcification_relevance_keeps_glycosyltransferase_plausible_with_housekeeping_noise(self):
        row = {
            "function_call": "GT47 exostosin-like membrane glycosyltransferase",
            "top_pfam_description": "Exostosin GT47 domain",
            "top_blast_description": "Putative PWWP domain-containing DNA repair factor 4",
            "top_pfam_summary": "Exostosin GT47 domain (2)",
            "top_blast_summary": "Putative PWWP domain-containing DNA repair factor 4 (1); Probable xyloglucan galactosyltransferase GT20 (1)",
            "signalp_fraction": 0.333,
            "tm_fraction": 1.0,
        }
        relevance, _ = infer_calcification_relevance(row)
        self.assertEqual(relevance, "plausible indirect")

    def test_infer_function_prefers_pfam_glycosyltransferase_over_unrelated_sparse_blast(self):
        row = {
            "top_pfam_description": "Exostosin GT47 domain",
            "top_blast_description": "Putative PWWP domain-containing DNA repair factor 4",
            "top_pfam_count": 2,
            "pfam_annotated_queries": 2,
            "top_blast_count": 1,
            "blast_annotated_queries": 3,
            "signalp_fraction": 0.333,
            "tm_fraction": 1.0,
            "glob_fraction": 0.0,
            "top_pfam_summary": "Exostosin GT47 domain (2)",
            "top_blast_summary": "Putative PWWP domain-containing DNA repair factor 4 (1); Probable xyloglucan galactosyltransferase GT20 (1)",
        }
        self.assertEqual(infer_function(row), "GT47 exostosin-like membrane glycosyltransferase")

    def test_infer_function_prefers_strong_full_length_helicase_hit_over_domain_only_pfam(self):
        row = {
            "top_pfam_description": "ATP-dependent RecD-like DNA helicase SH3 domain",
            "top_blast_description": "ATP-dependent RecD2 DNA helicase",
            "top_pfam_count": 10,
            "pfam_annotated_queries": 11,
            "top_blast_count": 11,
            "blast_annotated_queries": 11,
            "signalp_fraction": 0.091,
            "tm_fraction": 0.0,
            "glob_fraction": 1.0,
            "top_pfam_summary": "ATP-dependent RecD-like DNA helicase SH3 domain (10); UvrD-like helicase C-terminal domain (10)",
            "top_blast_summary": "ATP-dependent RecD2 DNA helicase (11)",
        }
        self.assertEqual(infer_function(row), "Likely ATP-dependent RecD2 DNA helicase")


if __name__ == "__main__":
    unittest.main()
