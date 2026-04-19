# Orthogroup annotation — calcifying / "loose" set

## Inputs used

- `*.fa` – 73 orthogroup FASTA files (per-OG protein sequences)
- `blast_out/OG*.blast.tsv` – BLAST vs UniProtKB/Swiss-Prot
  - columns: query, subject, subject-description, %identity, aln-length, e-value, bitscore
- `hmm_out/OG*.pfam.tbl` – HMMER `hmmscan` tabular output against Pfam
- `deeptmhmm_results/*.predicted_topologies.3line` – DeepTMHMM topology type per protein (`GLOB`, `SP`, `TM`, `SP+TM`, `BETA`)
- `signalp_results/signalp_batch_*.output.json` – SignalP 6.0 per-protein prediction (`Other` vs `Signal Peptide (Sec/SPI)` etc.)

## Generated files

- `annotate_ogs.py` – parser that aggregates the four evidence sources per OG
- `annotation_summary.tsv` – raw per-OG evidence
  - columns: `OG`, `n_sequences`, `blast_n_with_hit`, `top_blast_descriptions`, `pfam_n_with_domain`, `top_pfam_domains`, `signalp_counts`, `deeptmhmm_counts`
- `annotation_calcification.tsv` – interpreted per-OG annotation
  - columns: `OG`, `n_seq`, `primary_annotation`, `key_domains_or_hits`, `secretion_profile`, `calcification_relevance`, `rationale`

## Scoring rubric used for `calcification_relevance`

| Bucket | What it means | Count |
| --- | --- | --- |
| **HIGH** | Enzyme/protein family with documented role in biomineral organic matrix (sulfation, sulfated-polysaccharide / heparan-sulfate biosynthesis, matrix-modifying sulfatases). | 3 |
| **MEDIUM** | Secreted proteases, integrin-like adhesion proteins, matrix glycosyltransferases, cell-surface anchors, BMP/skeletogenic-signalling-adjacent kinases — plausible but indirect roles in mineral deposition. | 6 |
| **LOW-MEDIUM** | Enzyme chemistry (e.g. 2-OG/Fe(II) oxygenases, PP2C Ca2+/Mg2+ phosphatases, matrix-glycan-decorating transferases) that could touch biomineralization, but top BLAST hits don't support it directly. | 9 |
| **LOW** | Housekeeping / cytosolic / unrelated functions (DNA helicases, translation factors, chaperones, ubiquitin-system, secondary metabolism, pentapeptide-repeat proteins, etc.). | 33 |
| **UNKNOWN** | No assignable BLAST or Pfam evidence ("dark-matter" OGs). Some are notable because of strong SP/TM signal — flagged in the rationale. | 22 |

## Candidate calcification-relevant orthogroups

**HIGH relevance**

- `OG0017138` – Arylsulfatase B (Pfam `Sulfatase`, 14/16). Sulfatases act on sulfated polysaccharides abundant in biomineral organic matrices (coccolith polysaccharides, mollusc shell GAGs, echinoderm skeletal matrices).
- `OG0020703` – Heparan-sulfate 6-O-sulfotransferase (Pfam `Sulfotransfer_1/3`; BLAST HS6ST1). Strongly secreted. Direct role in building sulfated acidic polysaccharides that regulate CaCO3 nucleation.
- `OG0011061` – Exostosin-family GT47 glycosyltransferase. Type-II Golgi membrane topology (TM/SP+TM) consistent with heparan-sulfate / sulfated polysaccharide chain extension.

**MEDIUM relevance**

- `OG0009246` – Secreted collagen-VI / adhesive-plaque-matrix-like ECM candidates.
- `OG0009301` – PRIMA1-like proline-rich membrane anchor; SP+TM heavy. Plausible cell-surface tether.
- `OG0018986` – FG-GAP integrin-alpha-like propellers; cell-matrix divalent-cation-dependent adhesion.
- `OG0021347` – Secreted trypsin-family serine proteases; all members SP-positive. Matrix remodeling.
- `OG0023594` – GT8 glycosyltransferase (galacturonosyltransferase-like); builds acidic polysaccharide backbones.
- `OG0017305` – Stf0-like sulphotransferase; half the members secreted; potentially matrix-GAG sulfation.

**LOW-MEDIUM (watch-list)**

- `OG0023496` – BMP-2-inducible kinase (BMP2K/BIKE) and cyclin-G-associated kinase; tangential to BMP bone-signalling.
- `OG0021523` – PP2C Ser/Thr phosphatases (Ca2+/Mg2+-dependent); peripheral Ca2+-signalling link.
- `OG0011197` / `OG0026865` – 2-OG/Fe(II) dioxygenases and JmjC hydroxylases — same chemistry as collagen prolyl/lysyl hydroxylases but top hits don't point there.
- `OG0001332` – Hedgehog/Hint-domain morphogen-like (possible skeletogenic cell-fate role).
- `OG0009816`, `OG0024846` – mannosyl/hexosyltransferases; matrix-glycoprotein decoration.
- `OG0014155` – secreted nucleoside hydrolases (secretion is suggestive but chemistry isn't mineralization-specific).
- `OG0016203` – TM-anchored proteins with one collagen-IV hit.

**UNKNOWN but secretion-enriched — worth follow-up**

- `OG0020657` (SP 7/11), `OG0022355` (SP+TM 6/9), `OG0022473` (SP+TM 3/9), `OG0025059` (SP 5/7), `OG0010991` (small SP fraction), `OG0018519` (all TM), `OG0025015` (all TM).
These OGs lack BLAST/Pfam hits but show a clear secreted/membrane topology — the pattern expected of lineage-specific calcification matrix proteins. They are the most promising "dark-matter" candidates for structural modelling (e.g. ColabFold / ESMFold) or motif search (acidic / Ca2+-binding / EF-hand, GXGXG, pro-rich, DEAD-like acidic repeats).

## How to reproduce

```bash
cd /Users/xiaoyu/workspace/Orthogroups.calcifying_loose_fastas
python3 annotate_ogs.py    # regenerates annotation_summary.tsv
```

`annotation_calcification.tsv` is a curated interpretation layered on top of that raw evidence and is edited by hand.
