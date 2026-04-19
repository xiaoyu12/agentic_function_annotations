# Orthogroup Annotation & Calcification Relevance Report

Annotation of **71 orthogroups** (1705 sequences) from calcifying-haptophyte-enriched OGs.
Sources combined: BLAST vs SwissProt (E<1e-3 aware), HMMER-Pfam (E<1e-3), SignalP 6.0 secretion calls,
DeepTMHMM transmembrane predictions. Intermediate tables: `og_annotation_summary.tsv` / `og_details.md`.

**Species coverage** (JGI prefixes): Callep1130_1 (*Calcidiscus leptoporus*), Emihu1 (*Emiliania huxleyi*),
Gepoce1 (*Gephyrocapsa oceanica*), Hymcor1016_1, Ochro3194_1 (*Ochrosphaera*), Pleelo874_1 (*Pleurochrysis*),
Umbfol2878_1 (*Umbilicosphaera*), Chrlea393_1 / Chrsim136_1 (*Chrysotila*), plus a few non-calcifying
haptophytes (Isochrysis, Pavlova, Phaeocystis, Imantonia, Platychrysis) used as comparators.

## Calcification-relevance rubric

| Tier | Criterion |
|------|-----------|
| **HIGH** | Direct homolog to a known biomineralization / shell-matrix / coccolith protein, or strong "acidic secreted glycoprotein" signature (SP+, sulfation/acidic domain, high copy in calcifiers) |
| **MEDIUM** | Functional class implicated in biomineralization ECM (GT, sulfotransferase, sulfatase, GAG enzyme, integrin-like, collagen, mannosyltransferase, sulfate/pH signalling) or secreted + TM with no direct ortholog but topology fits |
| **LOW** | Plausible peripheral role (signalling kinases, proteases with SP, ubiquitin routing of matrix proteins) — candidates only |
| **UNLIKELY** | Housekeeping / metabolic with no link to calcification (DNA repair, translation, plastid catabolism, bacterial transferred genes) |
| **CANDIDATE (orphan)** | No BLAST/Pfam hit but strong secretory/TM topology conserved across calcifiers — novel candidate worth prioritizing |

---

## Master annotation table

| OG | n | SigP+ | TM+ | mean/max TMR | Best annotation | Putative function | Relevance |
|----|---|-------|-----|------|-----------------|-------------------|-----------|
| OG0000049 | 502 | 2 | 3 | 0.02/6 | FKBP15 / CCDC39 coiled-coil | Cytoskeletal coiled-coil / FK506-binding scaffold; weak hits only | UNLIKELY |
| OG0001332 | 91 | 1 | 76 | 0.90/2 | Warthog/Hedgehog (Hint module, 85/91) | Hint-domain autoprocessing proteins, Hh-like secreted signal (but here predominantly TM, SP−) | LOW (signalling) |
| OG0001976 | 73 | 4 | 3 | 0.07/2 | Pentapeptide-repeat (slr1819-like) | Pentapeptide-repeat proteins, cyanobacterial-like (likely HGT); no calcification link known | UNLIKELY |
| OG0002887 | 59 | 1 | 0 | 0/0 | PIF1-family ATP-dep DNA helicase | DNA replication/repair helicase | UNLIKELY |
| OG0003717 | 51 | 0 | 4 | 0.12/2 | — (no hits) | Orphan; weak TM signal | CANDIDATE (orphan, weak) |
| **OG0009246** | 29 | **9** | 4 | 0.14/1 | **Mytilus FP1 adhesive plaque / Col6α5** | **Secreted adhesive/matrix-like protein**; best hit to mussel foot protein FP1, plus collagen-VI; 31% SP+ | **HIGH** |
| OG0009301 | 29 | 21 | 15 | 0.55/2 | SHOCT (short C-term domain), PRIMA1 | Membrane-anchored small secreted proteins; ~72% SP+ and half TM — mostly orphan outside SHOCT | MEDIUM (secreted membrane) |
| OG0009816 | 28 | 2 | 6 | 0.36/5 | α-1,3-mannosyltransferase MNN14 | Mannosyl_trans3 family; N/O-glycan biosynthesis — could glycosylate coccolith polysaccharides | MEDIUM |
| OG0010177 | 28 | 0 | 0 | 0/0 | Ankyrin-repeat protein (Rickettsia RBE_0220-like) | Cytosolic ANK-repeat scaffold; possibly HGT | UNLIKELY |
| OG0010867 | 27 | 1 | 0 | 0/0 | Pentapeptide-repeat / SpkB-like | Pentapeptide-repeat (slr1819-type); no calcification link | UNLIKELY |
| OG0010991 | 27 | 4 | 1 | 0.07/2 | — (no hits) | Orphan with a few SP+ | CANDIDATE (weak) |
| **OG0011061** | 27 | **9** | **27** | **9.11/12** | Exostosin GT47 | **Exostosin-family glycosyltransferase (HS / polysaccharide chain-extension)**; every member TM; secretory pathway membrane enzyme — extends sulfated polysaccharides | **HIGH** |
| OG0011197 | 27 | 1 | 2 | 0.07/1 | 2OG-Fe(II) oxygenase (LFS-like) | Dioxygenase; secondary metabolism (plant-like) | UNLIKELY |
| OG0013037 | 24 | 5 | 6 | 0.25/1 | Hydrolase_4 (S33 serine aminopeptidase) | α/β-hydrolase; weak matches to giant muscle proteins are spurious | LOW |
| OG0014155 | 22 | 11 | 2 | 0.09/1 | Nucleoside hydrolase (NSH4 / RihA) | Purine/pyrimidine salvage; half SP+ — secreted nucleoside hydrolase (possible periplasmic) | UNLIKELY (metabolism) |
| OG0014250 | 22 | 15 | 1 | 0.05/1 | gatA — Glu-tRNA(Gln) amidotransferase A | Mitochondrial/organellar tRNA editing amidase | UNLIKELY |
| OG0015153 | 20 | 1 | 17 | 6.25/9 | piccolo (weak) | Multi-TM membrane protein, no clear ortholog; orphan TM family | CANDIDATE (multi-TM) |
| **OG0016203** | 18 | 1 | 15 | 6.39/12 | **Collagen α-2(IV) / Hexamerin 110** | **Multi-TM + collagen triple-helix-repeat proteins**; integral membrane collagen-like | **MEDIUM-HIGH** |
| **OG0016211** | 18 | 2 | **18** | 3.78/4 | SUR7 / PalI family | **4-TM membrane organizer** (membrane microdomain / MCC/eisosome in fungi); all members TM; plausible coccolith-vesicle membrane scaffold | **MEDIUM** |
| **OG0017138** | 16 | 2 | 0 | 0/0 | **Arylsulfatase B (Arsb)** | **Sulfatase; desulfation of sulfated GAGs** — directly implicated in coccolith-polysaccharide sulfation cycling; E=3.6e-95 to rat ARSB | **HIGH** |
| OG0017238 | 16 | 0 | 0 | 0/0 | tRNA ligase 1 (RNL) | tRNA splicing ligase | UNLIKELY |
| **OG0017305** | 16 | 1 | 0 | 0/0 | **Stf0 sulphotransferase** | **Carbohydrate/small-molecule sulfotransferase**; relevant to sulfated coccolith polysaccharide (CAP) biosynthesis | **MEDIUM** |
| OG0017361 | 16 | 1 | 2 | 0.12/1 | DJ-1/PfpI — Isonitrile hydratase | DJ-1/PfpI family cysteine hydrolase; stress response | UNLIKELY |
| OG0017362 | 16 | 7 | 2 | 0.12/1 | Red-chlorophyll-catabolite reductase | Plastid chlorophyll catabolism | UNLIKELY |
| OG0017884 | 15 | 0 | 4 | 1.27/8 | GNAT acetyltransferase (mshD-like) | Acyl-CoA acetyltransferase | UNLIKELY |
| OG0017903 | 15 | 0 | 0 | 0/0 | Flavohemoprotein (hmp, E=2.2e-106) | NO-detoxifying globin/FAD-oxidoreductase | UNLIKELY |
| OG0017965 | 15 | 0 | 0 | 0/0 | Kelch-repeat / NANM | Kelch β-propeller; possibly BTB-Kelch substrate adaptor | UNLIKELY |
| OG0018519 | 14 | 0 | 13 | 3.64/4 | — (no hits) | **Conserved multi-TM orphan** across calcifiers (13/14 TM+) | CANDIDATE (multi-TM) |
| OG0018986 | 13 | 7 | 10 | 1.92/8 | FG-GAP repeat / integrin-α-like / DEX1 | **Integrin-α-like FG-GAP β-propeller** with TM anchor; extracellular Ca²⁺/matrix-binding module (plant DEX1 is exine-deposition) | **MEDIUM-HIGH** |
| OG0019067 | 13 | 0 | 5 | 0.77/6 | LRR / NLRC3-like | Innate-immunity LRR / NLR-type receptor | UNLIKELY |
| OG0019078 | 13 | 0 | 0 | 0/0 | eEF1-γ (EEF1G; GST_N/C) | Translation elongation factor | UNLIKELY |
| OG0019174 | 13 | 5 | 0 | 0/0 | gatA amidase (mitochondrial paralog) | tRNA-Gln amidation | UNLIKELY |
| OG0019180 | 13 | 0 | 0 | 0/0 | — (no hits) | Orphan, no features | UNLIKELY |
| OG0019790 | 12 | 0 | 0 | 0/0 | HECT E3 ubiquitin ligase (hecd-1) + ANK | Protein-ubiquitination scaffold | LOW |
| OG0019817 | 12 | 8 | 1 | 0.08/1 | FAD/NADPH monooxygenase (auaG-like) | Secreted (67% SP+) flavin monooxygenase — secondary-metabolism (aurachin/phenoxazinone-like) | LOW |
| OG0019873 | 12 | 0 | 0 | 0/0 | ACC synthase (ACS11) | Aminotransferase I/II; plant ethylene biosynthesis, here likely cytosolic aminotransferase | UNLIKELY |
| OG0019887 | 12 | 0 | 0 | 0/0 | Intron maturase / MMR_HSR1 GTPase | Organellar intron splicing / ribosome-binding GTPase | UNLIKELY |
| OG0019917 | 12 | 1 | 3 | 0.33/2 | PKA regulatory subunit (pkaR) | cAMP-binding kinase regulatory subunit; **signalling (cAMP→PKA) — may regulate calcification** | LOW-MEDIUM |
| OG0020622 | 11 | 1 | 0 | 0/0 | RecD2 DNA helicase | Bacterial DNA repair helicase | UNLIKELY |
| OG0020657 | 11 | 3 | 0 | 0/0 | — (no hits) | Orphan, partially secreted | CANDIDATE (weak) |
| OG0020696 | 11 | 1 | 3 | 0.27/1 | — (no hits) | Orphan | CANDIDATE (weak) |
| OG0020700 | 11 | 0 | 0 | 0/0 | PDZ domain protein | PDZ-scaffold | LOW |
| **OG0020703** | 11 | **4** | 0 | 0/0 | **Heparan-sulfate 6-O-sulfotransferase 1 (HS6ST1)** | **Sulfotransferase that O-sulfates GAG sugars** — highly relevant to sulfated coccolith polysaccharide (CAP) maturation | **HIGH** |
| OG0020706 | 11 | 0 | 5 | 2.64/8 | — (no hits) | Orphan multi-TM | CANDIDATE (multi-TM) |
| OG0020728 | 11 | 0 | 0 | 0/0 | — (no hits) | Orphan, no features | UNLIKELY |
| OG0020729 | 11 | 0 | 0 | 0/0 | — (no hits) | Orphan, no features | UNLIKELY |
| OG0021347 | 10 | 5 | 0 | 0/0 | Trypsin-family serine protease | Secreted serine protease (50% SP+) — could mature matrix proteins | LOW-MEDIUM |
| OG0021406 | 10 | 0 | 0 | 0/0 | 2-dehydro-3-deoxy-D-gluconate 5-DH (kduD) | Short-chain dehydrogenase (bacterial-type sugar metabolism) | UNLIKELY |
| OG0021517 | 10 | 0 | 0 | 0/0 | F-box protein | SCF-complex substrate adaptor | LOW |
| OG0021520 | 10 | 0 | 1 | 1.40/14 | — (no hits) | Single highly-TM orphan; mostly no features | CANDIDATE (weak) |
| OG0021523 | 10 | 0 | 0 | 0/0 | PP2C protein phosphatase | Ser/Thr phosphatase; stress/Ca²⁺-signalling (plant-type) | LOW-MEDIUM |
| OG0021531 | 10 | 0 | 0 | 0/0 | DegP / HtrA serine protease + PDZ | Periplasmic/organellar quality-control protease | LOW |
| OG0022355 | 9 | 6 | 6 | 0.67/1 | DUF2237 (bacterial-conserved) | Secreted (67% SP+) and TM-anchored DUF protein | CANDIDATE (secreted) |
| OG0022381 | 9 | 0 | 4 | 0.44/1 | NADPH-cytochrome P450 reductase (CPR) | Microsomal electron donor for CYPs (ER-anchored) | UNLIKELY |
| OG0022455 | 9 | 0 | 0 | 0/0 | — (no hits) | Orphan, no features | UNLIKELY |
| OG0022473 | 9 | 3 | 3 | 0.33/1 | — (no hits) | Orphan, partial SP+/TM+ | CANDIDATE (weak) |
| OG0022474 | 9 | 0 | 0 | 0/0 | F-box protein | SCF substrate adaptor | LOW |
| OG0022492 | 9 | 0 | 0 | 0/0 | Myosin heavy chain kinase B / WD40 | Ser/Thr kinase + β-propeller | UNLIKELY |
| OG0022500 | 9 | 0 | 0 | 0/0 | — (no hits) | Orphan, no features | UNLIKELY |
| OG0022520 | 9 | 0 | 0 | 0/0 | — (no hits) | Orphan, no features | UNLIKELY |
| OG0022524 | 9 | 0 | 0 | 0/0 | Pentapeptide / MfpA-like | Pentapeptide-repeat (slr1819/MfpA) | UNLIKELY |
| OG0022528 | 9 | 0 | 0 | 0/0 | — (no hits) | Orphan | UNLIKELY |
| OG0023496 | 8 | 0 | 0 | 0/0 | BMP2K / Cyclin-G-associated kinase | Ser/Thr kinase (GAK/BMP2K-like; clathrin cycle) | LOW |
| OG0023566 | 8 | 3 | 1 | 0.12/1 | SEY1 / RHD3 GTPase | ER-membrane fusion GTPase | LOW |
| OG0023594 | 8 | 1 | 1 | 0.25/2 | Glycosyltransferase family 8 / GT-like | **GT8 glycosyltransferase**; candidate ECM/coccolith polysaccharide enzyme | **MEDIUM** |
| OG0023657 | 8 | 0 | 0 | 0/0 | — (no hits) | Orphan, no features | UNLIKELY |
| OG0023790 | 8 | 1 | 0 | 0/0 | Red-chlorophyll-catabolite reductase | Plastid chlorophyll catabolism (fragmentary, RCC_reductase) | UNLIKELY |
| OG0024846 | 7 | 1 | 0 | 0/0 | Hexosyltransferase MUCI70 / TOD1 | **Plant-type GT-like (TOD1/MUCI70)** — mucilage/cell-wall polysaccharide modification | **MEDIUM** |
| OG0024979 | 7 | 0 | 0 | 0/0 | — (no hits) | Orphan, no features | UNLIKELY |
| OG0025015 | 7 | 0 | **7** | 3.14/4 | — (no hits) | **Pan-calcifier multi-TM orphan** (7/7 TM+) | **CANDIDATE (strong orphan)** |
| OG0025016 | 7 | 0 | 0 | 0/0 | — (no hits) | Orphan, no features | UNLIKELY |
| OG0025059 | 7 | 5 | 0 | 0/0 | — (no hits) | **Pan-calcifier secreted orphan** (5/7 SP+) | **CANDIDATE (strong orphan)** |
| OG0026865 | 6 | 0 | 0 | 0/0 | JmjC / PSR / JMJD6 | JmjC-domain hydroxylase / lysyl-hydroxylase — could hydroxylate collagen-like substrates | LOW-MEDIUM |

---

## Top calcification-relevant candidates (by tier)

### HIGH relevance — matrix / sulfation / glycan-processing enzymes

1. **OG0009246 — Adhesive plaque / collagen-VI homolog** (29 members, 9 SP+).
   Best BLAST: *Mytilus coruscus* foot protein FP1 (E=2.1e-08) and mouse Col6α5 (E=2.9e-04).
   Biomineralizing bivalves use FP1-family proteins as acidic, DOPA-rich surface adhesives that
   nucleate CaCO₃. Strong signal-peptide enrichment is consistent with a secreted coccolith-
   matrix / adhesion protein. **Priority follow-up.**

2. **OG0011061 — Exostosin GT47 glycosyltransferase** (27 members, 9 SP+, 27/27 TM, mean 9 TMRs).
   Exostosin GT47 proteins polymerize glucuronic-acid / glucosamine backbones of heparan sulfate
   and xyloglucan. In coccolithophores, sulfated-coccolith-polysaccharides (CAPs) form the
   template for calcite nucleation — a GT47-family enzyme is a very strong candidate for
   CAP backbone synthesis. Every member is a polytopic secretory-pathway membrane protein.

3. **OG0017138 — Arylsulfatase B (ARSB)** (16 members, highly conserved E=3.6e-95).
   N-acetylgalactosamine-4-sulfatase; cleaves 4-sulfate from sulfated GAGs. In *E. huxleyi*,
   coccolith-associated polysaccharides carry uronic acids and sulfate esters — ARSB-family
   sulfatases likely remodel CAP sulfation state in the Golgi/coccolith-vesicle system.

4. **OG0020703 — Heparan-sulfate 6-O-sulfotransferase 1 (HS6ST1)** (11 members, 4 SP+).
   PAPS-dependent sulfotransferase acting on GAG sugar 6-OH. Direct candidate for installing
   sulfate groups on coccolith polysaccharides during their biogenesis.

### MEDIUM relevance — plausible ECM / membrane biomineralization components

- **OG0016203 — Membrane collagen-like / hexamerin** (15/18 TM, collagen triple-helix repeat).
  Polytopic ER/vacuolar membrane protein with extracellular collagen repeats — plausible
  transmembrane ECM scaffold around the coccolith vesicle (CV).
- **OG0016211 — SUR7 / PalI membrane microdomain** (18/18 TM). 4-TM organisers; in fungi SUR7
  defines plasma-membrane MCC eisosomes — a parallel role organising CV membrane microdomains
  is reasonable.
- **OG0018986 — Integrin-α FG-GAP-repeat + TM** (FG-GAP β-propeller, 10/13 TM, 7 SP+).
  FG-GAP repeats form a β-propeller that coordinates divalent cations (Mg²⁺/Ca²⁺/Mn²⁺) at the
  MIDAS site — a plausible Ca²⁺-binding extracellular adhesion / mineralization receptor.
- **OG0017305 — Stf0 sulphotransferase** (PAPS-dependent). Small-molecule/carbohydrate
  sulfotransferase; could participate in sulfated polysaccharide biosynthesis.
- **OG0009816 — α-1,3-mannosyltransferase (MNN14)**. Candidate for N-/O-glycan chain
  assembly on secreted coccolith-matrix glycoproteins.
- **OG0023594 — GT8 glycosyltransferase**. GT8 family builds α-galacturonan, pectin-like polymers,
  LPS O-chains — candidate for CAP / matrix polysaccharide synthesis.
- **OG0024846 — TOD1/MUCI70 hexosyltransferase-like**. Plant homologs are mucilage / cell-wall
  biosynthesis — candidate for external matrix polysaccharide assembly.
- **OG0009301 — SHOCT / PRIMA1-like short anchored proteins** (21 SP+, 15 TM). Small secreted,
  membrane-anchored proteins with no specific match beyond SHOCT — classic topology of
  lineage-specific biomineralization matrix proteins.
- **OG0019917 — PKA regulatory subunit (pkaR)**. cAMP-signalling regulator; cAMP/PKA control
  cellular pH and ion transport — may regulate calcification.
- **OG0021347 — Secreted trypsin-family protease** (5/10 SP+). Possibly processes matrix
  precursor proteins (analogous to astacin/BMP-1 in animal shell formation).
- **OG0026865 — JmjC PSR / lysyl-hydroxylase**. JmjC hydroxylases modify collagen-like
  substrates (P4H/LH analogs), potentially hydroxylating collagen-like OG0016203 products.

### CANDIDATE (orphan) — worth prioritising for experimental validation

OGs with **no BLAST/Pfam hits** but strong SP/TM conservation across calcifying species are
prime candidates for lineage-specific biomineralization proteins (the classic
"unknown secreted protein" pattern seen for GPA in *E. huxleyi*):

| OG | n | SigP+ | TM+ | Why interesting |
|----|---|-------|-----|-----------------|
| **OG0025015** | 7 | 0 | **7/7** | 100% multi-TM orphan present in 7 calcifying species |
| **OG0025059** | 7 | 5 | 0 | 5/7 secreted orphan across calcifiers |
| **OG0018519** | 14 | 0 | 13/14 | Multi-TM orphan, broadly conserved |
| OG0020706 | 11 | 0 | 5 | Multi-TM orphan (mean 2.6 TMRs, max 8) |
| OG0021520 | 10 | 0 | 1 (14 TMRs!) | One highly-TM (14 TMRs) orphan |
| OG0022355 | 9 | 6 | 6 | DUF2237 with majority SP+ and TM+ — uncharacterised secreted/membrane protein |
| OG0010991 | 27 | 4 | 1 | Partly secreted orphan |
| OG0003717 | 51 | 0 | 4 | Large orphan group (51 members), weak TM |
| OG0020657/0696/22473 | — | 1–3 SP+, some TM+ | Smaller orphans with secretion/TM signatures |

### LOW / UNLIKELY

The majority of the remaining OGs (FKBP, CCDC, tRNA-ligase, gatA amidase, RecD2, PIF1, eEF1γ,
kduD, flavohemoprotein, chlorophyll-catabolite reductases, ACC-synthase, etc.) look like
housekeeping, plastid metabolism or prokaryote-derived (likely HGT) activities with no
documented biomineralization role. Repeated **pentapeptide-repeat / slr1819 / MfpA** hits
(OG0001976, OG0010867, OG0022524) almost certainly represent bacterial-origin pentapeptide-
repeat proteins retained in haptophytes; there is no known calcification link.

A notable special case: **OG0001332 (Warthog/Hedgehog Hint-module, 91 members)**. The Hint
autoprocessing module is the hallmark of Hedgehog morphogens, but in haptophytes this OG is
predominantly TM, mostly SP-negative, and lacks a clean secreted-ligand architecture. It is
more likely a lineage-specific Hint-family protein than a direct biomineralization player;
left in LOW (signalling/cell-communication) pending expression data.

---

## Recommendations / next steps

1. **Prioritise OG0009246, OG0011061, OG0017138, OG0020703** for expression analysis under
   calcifying vs non-calcifying conditions — these are strong, biochemically-grounded
   candidates.
2. **OG0016203 / OG0016211 / OG0018986** together form a plausible "coccolith-vesicle
   membrane module" (polytopic membrane proteins with collagen repeats, eisosome-like
   organiser, and Ca²⁺-binding β-propeller).
3. **OG0025015, OG0025059, OG0018519, OG0022355** are orphan candidates — run pan-haptophyte
   presence/absence to confirm restriction to calcifying species, then structural prediction
   (AlphaFold) and InterProScan against the full Pfam+CDD+PANTHER set.
4. Consider re-running BLAST against a marine-eukaryote/coccolithophore-transcript database
   (MMETSP, EukProt) for the "no hit" OGs — SwissProt is too sparse for haptophyte-specific
   proteins.
5. For glycan-processing candidates (OG0009816, OG0011061, OG0020703, OG0023594, OG0024846),
   cross-check with known *E. huxleyi* CAP-biosynthesis mutants / transcriptome data.

---

## Files in this annotation directory

- `aggregate.py` — Python script that regenerates all summaries from the raw files
- `og_annotation_summary.tsv` — one row per OG with per-OG counts (BLAST / Pfam / SigP / TM)
- `og_details.md` — per-OG detailed top-hit breakdown (raw material for the table above)
- `OG_annotation_report.md` — **this report**
