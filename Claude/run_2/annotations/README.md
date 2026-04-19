# Orthogroup annotation — calcifying-loose set

Pipeline: `annotate_ogs.py` at the repo root. Inputs: `blast_out/<OG>.blast.tsv` (SwissProt), `hmm_out/<OG>.pfam.tbl` (Pfam-A), `signalp_results/signalp_batch_*.output.json` (SignalP 6.0), `deeptmhmm_results/deeptmhmm_batch_*.*` (DeepTMHMM).

## Outputs

- `og_summary.tsv` — master table, one row per OG, sorted by calcification score.
- `og_summary.md` — same data bucketed into High / Medium / Low / Background.
- `per_og/<OG>.txt` — detailed per-orthogroup summary: top SwissProt hits with member counts, Pfam domains with member counts, SignalP/TMHMM fractions, amino-acid composition, calcification tags, biomineralization flags, score.

## Scoring

- `calc_score` combines two signals:
  1. Keyword hits on SwissProt descriptions (top-5 per member) and Pfam descriptions using a curated list of calcification-relevant terms (carbonic anhydrase, bicarbonate/anion transporters, Ca pumps/channels, EF-hand, annexin, coccolith, acidic-rich matrix proteins, sulfation, glycosyltransferases, etc.).
  2. Biomineralization "candidate flags" based on the secretome/membrane pattern: high SignalP fraction, high fraction of uncharacterized-but-secreted members, multi-pass transmembrane families with poor BLAST coverage, acidic-rich or cysteine-rich composition.

- Buckets: High ≥10, Medium 5–9, Low 2–4, Background 0–1.

## How to read the highlights

### High / medium confidence — functionally linked to calcification

- **OG0017138 — Arylsulfatase B (Sulfatase Pfam, 14/16 members)**. Sulfated polysaccharides make up the coccolith-associated polysaccharide (CAP) fraction of emiliania calcification. Arylsulfatases drive sulfate ester turnover on coccolith acidic polysaccharides and on sulfated extracellular matrix. Highest-scoring OG in the set.

- **OG0023594 — Glycosyltransferase family 8**. GT8 family glycosyltransferases build acidic extracellular polysaccharides (pectin-like, galacturonans). These polysaccharides are leading components of calcification matrices in coccolithophores and in molluscan/foram shells.

- **OG0024846 — MUCI70 hexosyltransferase (TOD1/MUCI70 Pfam)**. TOD1/MUCI70 GT-like domain family, implicated in plant mucilage/cell-wall polysaccharide biosynthesis — a candidate for templating acidic polysaccharide of the coccolith organic base plate.

- **OG0020703 — Heparan-sulfate 6-O-sulfotransferase** (novel-secreted flag). Sulfotransferase activity on sulfated glycans, consistent with roles in acidic polysaccharide/matrix production; 36% of members carry an N-terminal signal peptide.

- **OG0017305 — Stf0 sulphotransferase** (uncharacterized in SwissProt, Pfam Sulphotransf 3x). Bacterial Stf0-family sulfotransferases build sulfated trehalose lipids; in algae, related enzymes would sulfonate surface polysaccharides relevant to coccolith matrix chemistry.

- **OG0009301 — PRIMA1/SHOCT-containing (secreted, multi-pass)**. 72% of members secrete (SignalP+) and ~50% are membrane-embedded. The PRIMA1 Pfam hit (Proline-rich membrane anchor 1) plus high novelty (only 1/29 members have BLAST support) marks this as a good candidate for a coccolithophore-specific secreted/membrane family.

### Low confidence but interesting secretome/membrane candidates

- **OG0014250 — Amidase (Amidase Pfam, 17x)**, 68% signal-peptide positive, 12.2% acidic — a secreted amidase family; role in matrix remodelling plausible but speculative.
- **OG0014155 — Ribonucleoside hydrolase (IU_nuc_hydro)**, 50% secreted. Unusual for this enzyme class to be in secretome; flagged.
- **OG0021347 — Trypsin family** (Trypsin Pfam, 10x), 50% secreted. Proteolytic trimming of acidic matrix proteins is a known biomineralization accessory activity.
- **OG0019817 — FAD-binding monooxygenase**, 67% secreted — unexpected for a redox enzyme.
- **OG0018986 — FG-GAP repeats / TcdB-toxin mid-N** — integrin-like repeats with multi-pass TM and high secretion; FG-GAP proteins mediate Ca²⁺-dependent cell-substrate adhesion.
- **OG0025059, OG0022473, OG0022355** — small, uncharacterized, secretion-biased orthogroups; prime candidates for taxon-specific biomineralization proteins ("dark matter" secretome).
- **OG0015153, OG0016211, OG0018519, OG0025015, OG0011061** — multi-pass TM families, mostly no BLAST support. Novel membrane-protein families; `SUR7` domain in OG0016211 relates to membrane compartmentalization / lipid raft organisation. These are plausible transporters/channels, but further annotation (InterPro / structural homology) is needed to commit to a Ca²⁺/HCO₃⁻ role.

### Likely not calcification-related

Orthogroups in the "Background" bucket are dominated by housekeeping families (pentapeptide repeats, DNA helicases, ankyrin repeats, F-box proteins, kinases, tRNA/ribosome machinery, GST/oxidoreductases, protein phosphatases, HECT E3 ligases). No carbonic anhydrase, no Ca²⁺-ATPase, no SLC4/SLC26, no annexin, and no coccolith-GPA/CAP hits were detected in the entire set. If the input was intended to enrich for calcification-related proteins, this set is enriched for **extracellular matrix biosynthesis** (sulfation, glycosylation) rather than for canonical Ca²⁺/bicarbonate transport machinery.

## Caveats

- SwissProt top-hit keyword screening misses calcification proteins whose closest SwissProt match is a divergent homolog. Re-BLAST against TrEMBL or MMETSP+calcifying-genomes would be informative.
- "no annotation" top-hit rows with high SignalP% but low BLAST coverage are the best candidates for novel, lineage-specific biomineralization proteins; a structural-homology search (Foldseek against PDB + AlphaFoldDB) on representatives would be the natural next step.
- Calcification scoring is keyword-based and conservative; scores are a triage signal, not a ground truth.
