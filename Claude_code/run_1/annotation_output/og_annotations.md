# Orthogroup annotation & calcification relevance

Annotations generated from BLASTp (SwissProt), HMMER/Pfam (Pfam-A), DeepTMHMM topology, and SignalP 6.0 secretion predictions for the `Orthogroups.calcifying_loose_fastas` set.

## Dataset

- 73 orthogroups
- 1705 total sequences across 22 source proteomes

Source proteomes (sequence count) - all haptophytes/coccolithophores and related chromists used as a *calcifying* query set:

- `Pleelo874_1`: 466
- `Callep1130_1`: 344
- `Umbfol2878_1`: 174
- `Ochro3194_1`: 165
- `Gepoce1`: 157
- `Emihu1`: 121
- `Hymcor1016_1`: 115
- `Chrlea393_1`: 54
- `IsochDm2_1`: 30
- `Isogal1`: 22
- `Chrsim136_1`: 16
- `Phaant1`: 16
- `Phaglo1`: 10
- `Platy1217_1`: 3
- `Pavlova6257_1`: 2
- `Pavlova5622_1`: 2
- `Pavlova2550_1`: 2
- `Imarot4780_1`: 2
- `Imarot4477_3`: 1
- `Imarot704_1`: 1
- `Psesor5268_1`: 1
- `Imarot2298_1`: 1

- **High** calcification relevance: 1
- **Medium** relevance: 10
- **Low / no evidence**: 62

## Relevance tiers

* **high** - PFAM or BLAST match to a canonical calcification / Ca²⁺ / biomineralization marker (carbonic anhydrase, EF-hand, annexin, SLC4/26, Ca²⁺-ATPase, galaxin, VWA/VWC, collagen, cadherin, chitin-binding, etc.).
* **medium** - weaker or indirect evidence (peripheral ECM/adhesion term or single keyword match).
* **low** - no keyword/domain evidence for calcification (still may be lineage-specific secreted proteins - see SP/TM columns).

## Summary table

| OG | Tier | Function | SP% | TM% | Rationale |
|---|---|---|---|---|---|
| OG0016203 | high | Collagen (x1) / Collagen alpha-2(IV) chain (x19) | 5.6 | 83.3 | PFAM:Collagen (collagen triple helix); keyword:COL (ECM/scaffold) / topology: membrane-bound (83% TM) |
| OG0009246 | medium | tRNA-synt_2 (x1) / Adhesive plaque matrix protein (x3) | 31.0 | 13.8 | keyword:COL (ECM/scaffold); keyword:adhesive-plaque (byssal/adhesive matrix protein) / topology: partially secreted (31% SP) |
| OG0011061 | medium | Exostosin_GT47 (x2) / Putative PWWP domain-containing DNA repair factor 4 (x3) | 33.3 | 100.0 | keyword:piccolo (Ca2+-sensor/active-zone protein) / topology: partially secreted (33% SP), membrane-bound (100% TM) |
| OG0014155 | medium | IU_nuc_hydro (x17) / Pyrimidine-specific ribonucleoside hydrolase RihA (x17) | 50.0 | 9.1 | secreted protein (candidate extracellular matrix) / topology: predominantly secreted (50% SP) |
| OG0014250 | medium | Amidase (x17) / Glutamyl-tRNA(Gln) amidotransferase subunit A (x47) | 68.2 | 4.5 | secreted protein (candidate extracellular matrix) / topology: predominantly secreted (68% SP) |
| OG0015153 | medium | Protein piccolo (x20) | 5.0 | 85.0 | keyword:piccolo (Ca2+-sensor/active-zone protein) / topology: membrane-bound (85% TM) |
| OG0017138 | medium | Sulfatase (x14) / Arylsulfatase B (x44) | 12.5 | 0.0 | keyword:arylsulfatase (GAG/sulfated polysaccharide processing) |
| OG0019817 | medium | FAD_binding_3 (x9) / Zeaxanthin epoxidase, chloroplastic (x12) | 66.7 | 8.3 | secreted protein (candidate extracellular matrix) / topology: predominantly secreted (67% SP) |
| OG0020703 | medium | Sulfotransfer_3 (x1) / Heparan-sulfate 6-O-sulfotransferase 1 (x2) | 36.4 | 0.0 | keyword:sulfotransferase (sulfation of ECM/GAGs); keyword:heparan (HS proteoglycan) / topology: partially secreted (36% SP) |
| OG0021347 | medium | Trypsin (x10) / Trypsin-1 (x7) | 50.0 | 0.0 | secreted protein (candidate extracellular matrix) / topology: predominantly secreted (50% SP) |
| OG0025059 | medium | no detectable homology or domain | 71.4 | 0.0 | secreted protein (candidate extracellular matrix) / topology: predominantly secreted (71% SP) |
| OG0000049 | low | FKBP15 (x7) / FK506-binding protein 15 (x10) | 0.4 | 0.6 | no direct evidence |
| OG0001332 | low | Hint (x85) / Protein hedgehog (x81) | 1.1 | 83.5 | no direct evidence / topology: membrane-bound (84% TM) |
| OG0001976 | low | Pentapeptide (x72) / Uncharacterized protein in mobD 3'region (x220) | 5.5 | 4.1 | no direct evidence |
| OG0002887 | low | PIF1 (x18) / ATP-dependent DNA helicase PIF1 (x14) | 1.7 | 0.0 | no direct evidence |
| OG0003717 | low | no detectable homology or domain | 0.0 | 7.8 | no direct evidence |
| OG0009301 | low | SHOCT (x4) / DNA-directed RNA polymerase II subunit RPB1 (x3) | 72.4 | 51.7 | no direct evidence / topology: predominantly secreted (72% SP), membrane-bound (52% TM) |
| OG0009816 | low | Mannosyl_trans3 (x15) / Putative alpha-1,3-mannosyltransferase MNN14 (x10) | 7.1 | 21.4 | no direct evidence / topology: partial TM (21% TM) |
| OG0010177 | low | Ank_2 (x24) / Ankyrin repeat domain-containing protein 17 (x91) | 0.0 | 0.0 | no direct evidence |
| OG0010867 | low | Pentapeptide (x27) / Uncharacterized protein in mobD 3'region (x79) | 3.7 | 0.0 | no direct evidence |
| OG0010991 | low | no detectable homology or domain | 14.8 | 3.7 | no direct evidence |
| OG0011197 | low | 2OG-FeII_Oxy (x3) / RAP domain-containing protein, chloroplastic (x2) | 3.7 | 7.4 | no direct evidence |
| OG0013037 | low | Hydrolase_4 (x1) / Protein Obscurin (x1) | 20.8 | 25.0 | no direct evidence / topology: partially secreted (21% SP), partial TM (25% TM) |
| OG0016211 | low | SUR7 (x1) | 11.1 | 100.0 | no direct evidence / topology: membrane-bound (100% TM) |
| OG0017238 | low | tRNA ligase 1 (x11) | 0.0 | 0.0 | no direct evidence |
| OG0017305 | low | Sulphotransf (x3) | 6.2 | 0.0 | no direct evidence |
| OG0017361 | low | DJ-1_PfpI (x16) / Isonitrile hydratase (x18) | 6.2 | 12.5 | no direct evidence |
| OG0017362 | low | RCC_reductase (x3) / Red chlorophyll catabolite reductase 1, chloroplastic (x2) | 43.8 | 12.5 | no direct evidence / topology: partially secreted (44% SP) |
| OG0017884 | low | Acetyltransf_1 (x7) / Mycothiol acetyltransferase (x2) | 0.0 | 26.7 | no direct evidence / topology: partial TM (27% TM) |
| OG0017903 | low | FAD_binding_6 (x10) / Flavohemoprotein (x35) | 0.0 | 0.0 | no direct evidence |
| OG0017965 | low | Kelch_1 (x5) / N-acetylneuraminate epimerase (x4) | 0.0 | 0.0 | no direct evidence |
| OG0018519 | low | no detectable homology or domain | 0.0 | 92.9 | no direct evidence / topology: membrane-bound (93% TM) |
| OG0018986 | low | FG-GAP_3 (x12) / Metabotropic glutamate receptor-like protein E (x2) | 53.8 | 76.9 | no direct evidence / topology: predominantly secreted (54% SP), membrane-bound (77% TM) |
| OG0019067 | low | LRR_6 (x12) / NLR family CARD domain-containing protein 3 (x51) | 0.0 | 38.5 | no direct evidence / topology: partial TM (38% TM) |
| OG0019078 | low | GST_N (x11) / Elongation factor 1-gamma (x14) | 0.0 | 0.0 | no direct evidence |
| OG0019174 | low | Amidase (x13) / Glutamyl-tRNA(Gln) amidotransferase subunit A (x29) | 38.5 | 0.0 | no direct evidence / topology: partially secreted (38% SP) |
| OG0019180 | low | no detectable homology or domain | 0.0 | 0.0 | no direct evidence |
| OG0019790 | low | HECT (x9) / E3 ubiquitin-protein ligase TRIP12 (x9) | 0.0 | 0.0 | no direct evidence |
| OG0019873 | low | Aminotran_1_2 (x10) / 1-aminocyclopropane-1-carboxylate synthase 11 (x10) | 0.0 | 0.0 | no direct evidence |
| OG0019887 | low | MMR_HSR1 (x4) / Mitochondrial ribosome-associated GTPase 1 (x4) | 0.0 | 0.0 | no direct evidence |
| OG0019917 | low | cNMP_binding (x2) / cAMP-dependent protein kinase regulatory subunit (x8) | 8.3 | 25.0 | no direct evidence / topology: partial TM (25% TM) |
| OG0020622 | low | AAA_30 (x10) / ATP-dependent RecD2 DNA helicase (x26) | 9.1 | 0.0 | no direct evidence |
| OG0020657 | low | no detectable homology or domain | 27.3 | 0.0 | no direct evidence / topology: partially secreted (27% SP) |
| OG0020696 | low | no detectable homology or domain | 9.1 | 27.3 | no direct evidence / topology: partial TM (27% TM) |
| OG0020700 | low | PDZ_2 (x3) | 0.0 | 0.0 | no direct evidence |
| OG0020706 | low | no detectable homology or domain | 0.0 | 45.5 | no direct evidence / topology: partial TM (45% TM) |
| OG0020728 | low | no detectable homology or domain | 0.0 | 0.0 | no direct evidence |
| OG0020729 | low | no detectable homology or domain | 0.0 | 0.0 | no direct evidence |
| OG0021406 | low | adh_short (x10) / 2-dehydro-3-deoxy-D-gluconate 5-dehydrogenase (x17) | 0.0 | 0.0 | no direct evidence |
| OG0021517 | low | F-box-like (x8) / F-box only protein 3 (x3) | 0.0 | 0.0 | no direct evidence |
| OG0021520 | low | no detectable homology or domain | 0.0 | 10.0 | no direct evidence |
| OG0021523 | low | PP2C (x6) / Probable protein phosphatase 2C 64 (x5) | 0.0 | 0.0 | no direct evidence |
| OG0021531 | low | PDZ_6 (x10) / Probable periplasmic serine endoprotease DegP-like (x6) | 0.0 | 0.0 | no direct evidence |
| OG0022355 | low | DUF2237 (x7) | 66.7 | 66.7 | no direct evidence / topology: predominantly secreted (67% SP), membrane-bound (67% TM) |
| OG0022381 | low | FAD_binding_1 (x8) / NADPH--cytochrome P450 reductase (x15) | 0.0 | 44.4 | no direct evidence / topology: partial TM (44% TM) |
| OG0022455 | low | no detectable homology or domain | 0.0 | 0.0 | no direct evidence |
| OG0022473 | low | no detectable homology or domain | 33.3 | 33.3 | no direct evidence / topology: partially secreted (33% SP), partial TM (33% TM) |
| OG0022474 | low | F-box-like (x8) | 0.0 | 0.0 | no direct evidence |
| OG0022492 | low | Beta-prop_TEP1_2nd (x7) / Myosin heavy chain kinase B (x5) | 0.0 | 0.0 | no direct evidence |
| OG0022500 | low | no detectable homology or domain | 0.0 | 0.0 | no direct evidence |
| OG0022520 | low | no detectable homology or domain | 0.0 | 0.0 | no direct evidence |
| OG0022524 | low | Pentapeptide (x9) / Uncharacterized protein in mobD 3'region (x71) | 0.0 | 0.0 | no direct evidence |
| OG0022528 | low | no detectable homology or domain | 0.0 | 0.0 | no direct evidence |
| OG0023496 | low | Pkinase (x8) / Cyclin-G-associated kinase (x12) | 0.0 | 0.0 | no direct evidence |
| OG0023566 | low | RHD3_GTPase (x3) / Protein SEY1 (x5) | 37.5 | 12.5 | no direct evidence / topology: partially secreted (38% SP) |
| OG0023594 | low | Glyco_transf_8 (x6) / Ankyrin repeat, SAM and basic leucine zipper domain-containing protein 1 (x6) | 12.5 | 12.5 | no direct evidence |
| OG0023657 | low | no detectable homology or domain | 0.0 | 0.0 | no direct evidence |
| OG0023790 | low | RCC_reductase (x6) / Red chlorophyll catabolite reductase (Fragment) (x3) | 12.5 | 0.0 | no direct evidence |
| OG0024846 | low | TOD1_MUCI70 (x4) / Probable hexosyltransferase MUCI70 (x1) | 14.3 | 0.0 | no direct evidence |
| OG0024979 | low | no detectable homology or domain | 0.0 | 0.0 | no direct evidence |
| OG0025015 | low | no detectable homology or domain | 0.0 | 100.0 | no direct evidence / topology: membrane-bound (100% TM) |
| OG0025016 | low | no detectable homology or domain | 0.0 | 0.0 | no direct evidence |
| OG0026865 | low | JmjC (x5) / Bifunctional arginine demethylase and lysyl-hydroxylase JMJD6 (x7) | 0.0 | 0.0 | no direct evidence |

## Per-OG annotations (sorted by relevance)

### OG0016203  - _HIGH_

- sequences: 18; mean length: 441 aa
- signal peptide: 5.6%; ≥1 TM region: 83.3% (mean TMRs 6.39)
- top PFAM domains: Collagen(1); RCC1(1); WD40_RLD(1); RCC1_2(1)
- top BLAST subjects: Collagen alpha-2(IV) chain (x19); Collagen alpha-1(IV) chain (x18); Hexamerin 110 (x8)
- **function**: Collagen (x1) / Collagen alpha-2(IV) chain (x19)
- calcification rationale: PFAM:Collagen (collagen triple helix); keyword:COL (ECM/scaffold) | topology: membrane-bound (83% TM)

### OG0009246  - _MEDIUM_

- sequences: 29; mean length: 342 aa
- signal peptide: 31.0%; ≥1 TM region: 13.8% (mean TMRs 0.14)
- top PFAM domains: tRNA-synt_2(1)
- top BLAST subjects: Adhesive plaque matrix protein (x3); Collagen alpha-5(VI) chain (x1)
- **function**: tRNA-synt_2 (x1) / Adhesive plaque matrix protein (x3)
- calcification rationale: keyword:COL (ECM/scaffold); keyword:adhesive-plaque (byssal/adhesive matrix protein) | topology: partially secreted (31% SP)

### OG0011061  - _MEDIUM_

- sequences: 27; mean length: 607 aa
- signal peptide: 33.3%; ≥1 TM region: 100.0% (mean TMRs 9.11)
- top PFAM domains: Exostosin_GT47(2)
- top BLAST subjects: Putative PWWP domain-containing DNA repair factor 4 (x3); Protein piccolo (x3); Titin homolog (x1)
- **function**: Exostosin_GT47 (x2) / Putative PWWP domain-containing DNA repair factor 4 (x3)
- calcification rationale: keyword:piccolo (Ca2+-sensor/active-zone protein) | topology: partially secreted (33% SP), membrane-bound (100% TM)

### OG0014155  - _MEDIUM_

- sequences: 22; mean length: 364 aa
- signal peptide: 50.0%; ≥1 TM region: 9.1% (mean TMRs 0.09)
- top PFAM domains: IU_nuc_hydro(17)
- top BLAST subjects: Pyrimidine-specific ribonucleoside hydrolase RihA (x17); Nucleoside hydrolase 3 (x16); Nucleoside hydrolase 4 (x8)
- **function**: IU_nuc_hydro (x17) / Pyrimidine-specific ribonucleoside hydrolase RihA (x17)
- calcification rationale: secreted protein (candidate extracellular matrix) | topology: predominantly secreted (50% SP)

### OG0014250  - _MEDIUM_

- sequences: 22; mean length: 584 aa
- signal peptide: 68.2%; ≥1 TM region: 4.5% (mean TMRs 0.05)
- top PFAM domains: Amidase(17)
- top BLAST subjects: Glutamyl-tRNA(Gln) amidotransferase subunit A (x47); Putative amidase AmiD (x2); Glutamyl-tRNA(Gln) amidotransferase subunit A 1 (x2)
- **function**: Amidase (x17) / Glutamyl-tRNA(Gln) amidotransferase subunit A (x47)
- calcification rationale: secreted protein (candidate extracellular matrix) | topology: predominantly secreted (68% SP)

### OG0015153  - _MEDIUM_

- sequences: 20; mean length: 688 aa
- signal peptide: 5.0%; ≥1 TM region: 85.0% (mean TMRs 6.25)
- top PFAM domains: -
- top BLAST subjects: Protein piccolo (x20); Titin (x3); Subtelomeric hrmA-associated cluster protein AFUA_5G14930 (x2)
- **function**: Protein piccolo (x20)
- calcification rationale: keyword:piccolo (Ca2+-sensor/active-zone protein) | topology: membrane-bound (85% TM)

### OG0017138  - _MEDIUM_

- sequences: 16; mean length: 516 aa
- signal peptide: 12.5%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Sulfatase(14); Polysacc_lyase(1); Kelch_1(1); Kelch_2(1); Kelch_FKB95(1)
- top BLAST subjects: Arylsulfatase B (x44); Arylsulfatase I (x6); Arylsulfatase J (x1)
- **function**: Sulfatase (x14) / Arylsulfatase B (x44)
- calcification rationale: keyword:arylsulfatase (GAG/sulfated polysaccharide processing)

### OG0019817  - _MEDIUM_

- sequences: 12; mean length: 369 aa
- signal peptide: 66.7%; ≥1 TM region: 8.3% (mean TMRs 0.08)
- top PFAM domains: FAD_binding_3(9); NAD_binding_8(9)
- top BLAST subjects: Zeaxanthin epoxidase, chloroplastic (x12); Aurachin C monooxygenase/isomerase (x7); 3-amino-4-hydroxybenzoate 2-monooxygenase (x3)
- **function**: FAD_binding_3 (x9) / Zeaxanthin epoxidase, chloroplastic (x12)
- calcification rationale: secreted protein (candidate extracellular matrix) | topology: predominantly secreted (67% SP)

### OG0020703  - _MEDIUM_

- sequences: 11; mean length: 325 aa
- signal peptide: 36.4%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Sulfotransfer_3(1); Sulfotransfer_1(1)
- top BLAST subjects: Heparan-sulfate 6-O-sulfotransferase 1 (x2); Heparan-sulfate 6-O-sulfotransferase 2 (x2)
- **function**: Sulfotransfer_3 (x1) / Heparan-sulfate 6-O-sulfotransferase 1 (x2)
- calcification rationale: keyword:sulfotransferase (sulfation of ECM/GAGs); keyword:heparan (HS proteoglycan) | topology: partially secreted (36% SP)

### OG0021347  - _MEDIUM_

- sequences: 10; mean length: 313 aa
- signal peptide: 50.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Trypsin(10); Trypsin_2(3); DUF1986(1)
- top BLAST subjects: Trypsin-1 (x7); Trypsin-4 (x4); Trypsin 3A1 (x3)
- **function**: Trypsin (x10) / Trypsin-1 (x7)
- calcification rationale: secreted protein (candidate extracellular matrix) | topology: predominantly secreted (50% SP)

### OG0025059  - _MEDIUM_

- sequences: 7; mean length: 201 aa
- signal peptide: 71.4%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: secreted protein (candidate extracellular matrix) | topology: predominantly secreted (71% SP)

### OG0000049  - _LOW_

- sequences: 502; mean length: 673 aa
- signal peptide: 0.4%; ≥1 TM region: 0.6% (mean TMRs 0.02)
- top PFAM domains: FKBP15(7); CCDC39(2); ODF2_C(2); T4SS_LegC2C7(1); MAD(1)
- top BLAST subjects: FK506-binding protein 15 (x10); Reticulocyte-binding protein homolog 2a (x6); Myosin-10 (x3)
- **function**: FKBP15 (x7) / FK506-binding protein 15 (x10)
- calcification rationale: no direct evidence

### OG0001332  - _LOW_

- sequences: 91; mean length: 289 aa
- signal peptide: 1.1%; ≥1 TM region: 83.5% (mean TMRs 0.90)
- top PFAM domains: Hint(85); Trypsin(2); Hint_2(2); Recep_L_domain(1); Trypsin_2(1)
- top BLAST subjects: Protein hedgehog (x81); Warthog protein 1 (x76); Warthog protein 8 (x7)
- **function**: Hint (x85) / Protein hedgehog (x81)
- calcification rationale: no direct evidence | topology: membrane-bound (84% TM)

### OG0001976  - _LOW_

- sequences: 73; mean length: 302 aa
- signal peptide: 5.5%; ≥1 TM region: 4.1% (mean TMRs 0.07)
- top PFAM domains: Pentapeptide(72); Pentapeptide_4(70); Pentapeptide_3(70); Decapeptide(39); LD_SV2(24)
- top BLAST subjects: Uncharacterized protein in mobD 3'region (x220); Uncharacterized protein slr1819 (x203); Secreted effector protein PipB2 (x119)
- **function**: Pentapeptide (x72) / Uncharacterized protein in mobD 3'region (x220)
- calcification rationale: no direct evidence

### OG0002887  - _LOW_

- sequences: 59; mean length: 281 aa
- signal peptide: 1.7%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: PIF1(18); Pif1_2B_dom(17); AAA_30(7); DnaJ(3); Helitron_like_N(2)
- top BLAST subjects: ATP-dependent DNA helicase PIF1 (x14); ATP-dependent DNA helicase Pif1 (x11); ATP-dependent DNA helicase PIF4 (x4)
- **function**: PIF1 (x18) / ATP-dependent DNA helicase PIF1 (x14)
- calcification rationale: no direct evidence

### OG0003717  - _LOW_

- sequences: 51; mean length: 338 aa
- signal peptide: 0.0%; ≥1 TM region: 7.8% (mean TMRs 0.12)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0009301  - _LOW_

- sequences: 29; mean length: 219 aa
- signal peptide: 72.4%; ≥1 TM region: 51.7% (mean TMRs 0.55)
- top PFAM domains: SHOCT(4); PRIMA1(1)
- top BLAST subjects: DNA-directed RNA polymerase II subunit RPB1 (x3)
- **function**: SHOCT (x4) / DNA-directed RNA polymerase II subunit RPB1 (x3)
- calcification rationale: no direct evidence | topology: predominantly secreted (72% SP), membrane-bound (52% TM)

### OG0009816  - _LOW_

- sequences: 28; mean length: 502 aa
- signal peptide: 7.1%; ≥1 TM region: 21.4% (mean TMRs 0.36)
- top PFAM domains: Mannosyl_trans3(15)
- top BLAST subjects: Putative alpha-1,3-mannosyltransferase MNN14 (x10); Putative alpha-1,3-mannosyltransferase MNN12 (x7); Putative alpha-1,3-mannosyltransferase MNN13 (x7)
- **function**: Mannosyl_trans3 (x15) / Putative alpha-1,3-mannosyltransferase MNN14 (x10)
- calcification rationale: no direct evidence | topology: partial TM (21% TM)

### OG0010177  - _LOW_

- sequences: 28; mean length: 315 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Ank_2(24); Ank(24); Ank_5(24); Ank_4(24); Ank_3(24)
- top BLAST subjects: Ankyrin repeat domain-containing protein 17 (x91); Putative ankyrin repeat protein RBE_0220 (x48); Ankyrin repeat domain-containing protein 29 (x47)
- **function**: Ank_2 (x24) / Ankyrin repeat domain-containing protein 17 (x91)
- calcification rationale: no direct evidence

### OG0010867  - _LOW_

- sequences: 27; mean length: 340 aa
- signal peptide: 3.7%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Pentapeptide(27); Pentapeptide_3(26); Pentapeptide_4(25); Decapeptide(9); LD_SV2(5)
- top BLAST subjects: Uncharacterized protein in mobD 3'region (x79); Uncharacterized protein slr1819 (x70); Uncharacterized protein slr1152 (x42)
- **function**: Pentapeptide (x27) / Uncharacterized protein in mobD 3'region (x79)
- calcification rationale: no direct evidence

### OG0010991  - _LOW_

- sequences: 27; mean length: 337 aa
- signal peptide: 14.8%; ≥1 TM region: 3.7% (mean TMRs 0.07)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0011197  - _LOW_

- sequences: 27; mean length: 344 aa
- signal peptide: 3.7%; ≥1 TM region: 7.4% (mean TMRs 0.07)
- top PFAM domains: 2OG-FeII_Oxy(3)
- top BLAST subjects: RAP domain-containing protein, chloroplastic (x2); Kihadalactone A synthase LFS (x2); Gibberellin 2-beta-dioxygenase 7 (x1)
- **function**: 2OG-FeII_Oxy (x3) / RAP domain-containing protein, chloroplastic (x2)
- calcification rationale: no direct evidence

### OG0013037  - _LOW_

- sequences: 24; mean length: 451 aa
- signal peptide: 20.8%; ≥1 TM region: 25.0% (mean TMRs 0.25)
- top PFAM domains: Hydrolase_4(1)
- top BLAST subjects: Protein Obscurin (x1); Titin (x1)
- **function**: Hydrolase_4 (x1) / Protein Obscurin (x1)
- calcification rationale: no direct evidence | topology: partially secreted (21% SP), partial TM (25% TM)

### OG0016211  - _LOW_

- sequences: 18; mean length: 234 aa
- signal peptide: 11.1%; ≥1 TM region: 100.0% (mean TMRs 3.78)
- top PFAM domains: SUR7(1)
- top BLAST subjects: -
- **function**: SUR7 (x1)
- calcification rationale: no direct evidence | topology: membrane-bound (100% TM)

### OG0017238  - _LOW_

- sequences: 16; mean length: 552 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: tRNA ligase 1 (x11)
- **function**: tRNA ligase 1 (x11)
- calcification rationale: no direct evidence

### OG0017305  - _LOW_

- sequences: 16; mean length: 286 aa
- signal peptide: 6.2%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Sulphotransf(3)
- top BLAST subjects: -
- **function**: Sulphotransf (x3)
- calcification rationale: no direct evidence

### OG0017361  - _LOW_

- sequences: 16; mean length: 321 aa
- signal peptide: 6.2%; ≥1 TM region: 12.5% (mean TMRs 0.12)
- top PFAM domains: DJ-1_PfpI(16); HTH_18(1); HTH_AraC(1)
- top BLAST subjects: Isonitrile hydratase (x18); Isonitrile hydratase-like protein xanA (x6); HTH-type transcriptional regulator GbdR (x6)
- **function**: DJ-1_PfpI (x16) / Isonitrile hydratase (x18)
- calcification rationale: no direct evidence

### OG0017362  - _LOW_

- sequences: 16; mean length: 571 aa
- signal peptide: 43.8%; ≥1 TM region: 12.5% (mean TMRs 0.12)
- top PFAM domains: RCC_reductase(3); Exo_endo_phos(1); Abhydrolase_5(1)
- top BLAST subjects: Red chlorophyll catabolite reductase 1, chloroplastic (x2); DNA-(apurinic or apyrimidinic site) endonuclease (x2); DNA repair nuclease/redox regulator APEX1 (x1)
- **function**: RCC_reductase (x3) / Red chlorophyll catabolite reductase 1, chloroplastic (x2)
- calcification rationale: no direct evidence | topology: partially secreted (44% SP)

### OG0017884  - _LOW_

- sequences: 15; mean length: 384 aa
- signal peptide: 0.0%; ≥1 TM region: 26.7% (mean TMRs 1.27)
- top PFAM domains: Acetyltransf_1(7); Acetyltransf_7(7); RNF34L-like_3rd(2)
- top BLAST subjects: Mycothiol acetyltransferase (x2)
- **function**: Acetyltransf_1 (x7) / Mycothiol acetyltransferase (x2)
- calcification rationale: no direct evidence | topology: partial TM (27% TM)

### OG0017903  - _LOW_

- sequences: 15; mean length: 364 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: FAD_binding_6(10); Globin(5); NAD_binding_1(2); NPF(1); Protoglobin(1)
- top BLAST subjects: Flavohemoprotein (x35); Flavohemoprotein B (x5); Flavohemoprotein A (x2)
- **function**: FAD_binding_6 (x10) / Flavohemoprotein (x35)
- calcification rationale: no direct evidence

### OG0017965  - _LOW_

- sequences: 15; mean length: 326 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Kelch_1(5); Kelch_2(5); NANM(5); Kelch_KLHDC2_KLHL20_DRC7(5); Beta-prop_Calicin(5)
- top BLAST subjects: N-acetylneuraminate epimerase (x4); Kelch repeat-containing protein At3g27220 (x3); Kelch-like ECH-associated protein 1 (x3)
- **function**: Kelch_1 (x5) / N-acetylneuraminate epimerase (x4)
- calcification rationale: no direct evidence

### OG0018519  - _LOW_

- sequences: 14; mean length: 158 aa
- signal peptide: 0.0%; ≥1 TM region: 92.9% (mean TMRs 3.64)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence | topology: membrane-bound (93% TM)

### OG0018986  - _LOW_

- sequences: 13; mean length: 782 aa
- signal peptide: 53.8%; ≥1 TM region: 76.9% (mean TMRs 1.92)
- top PFAM domains: FG-GAP_3(12); FG-GAP(12); TcdB_toxin_midN(5); 7tm_3(2); BBS2_Mid(1)
- top BLAST subjects: Metabotropic glutamate receptor-like protein E (x2); Protein DEFECTIVE IN EXINE FORMATION 1 (x1); Gamma-aminobutyric acid type B receptor subunit 2 (x1)
- **function**: FG-GAP_3 (x12) / Metabotropic glutamate receptor-like protein E (x2)
- calcification rationale: no direct evidence | topology: predominantly secreted (54% SP), membrane-bound (77% TM)

### OG0019067  - _LOW_

- sequences: 13; mean length: 384 aa
- signal peptide: 0.0%; ≥1 TM region: 38.5% (mean TMRs 0.77)
- top PFAM domains: LRR_6(12); LRR_TMOD-LMOD(7); LRR_4(1)
- top BLAST subjects: NLR family CARD domain-containing protein 3 (x51); Protein NLRC3 (x39); Nucleotide-binding oligomerization domain-containing protein 1 (x15)
- **function**: LRR_6 (x12) / NLR family CARD domain-containing protein 3 (x51)
- calcification rationale: no direct evidence | topology: partial TM (38% TM)

### OG0019078  - _LOW_

- sequences: 13; mean length: 366 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: GST_N(11); GST_C(10); GST_N_2(10); GST_N_3(10); GST_C_2(9)
- top BLAST subjects: Elongation factor 1-gamma (x14); Elongation factor 1-gamma 2 (x8); Elongation factor 1-gamma 1 (x7)
- **function**: GST_N (x11) / Elongation factor 1-gamma (x14)
- calcification rationale: no direct evidence

### OG0019174  - _LOW_

- sequences: 13; mean length: 668 aa
- signal peptide: 38.5%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Amidase(13)
- top BLAST subjects: Glutamyl-tRNA(Gln) amidotransferase subunit A (x29); Glutamyl-tRNA(Gln) amidotransferase subunit A, mitochondrial (x9); Glutamyl-tRNA(Gln) amidotransferase subunit A 1 (x1)
- **function**: Amidase (x13) / Glutamyl-tRNA(Gln) amidotransferase subunit A (x29)
- calcification rationale: no direct evidence | topology: partially secreted (38% SP)

### OG0019180  - _LOW_

- sequences: 13; mean length: 160 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0019790  - _LOW_

- sequences: 12; mean length: 1103 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: HECT(9); Ank_2(4); Ank(4); Ank_3(4); Ank_4(4)
- top BLAST subjects: E3 ubiquitin-protein ligase TRIP12 (x9); E3 ubiquitin-protein ligase hecd-1 (x7); E3 ubiquitin-protein ligase UPL3 (x7)
- **function**: HECT (x9) / E3 ubiquitin-protein ligase TRIP12 (x9)
- calcification rationale: no direct evidence

### OG0019873  - _LOW_

- sequences: 12; mean length: 381 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Aminotran_1_2(10)
- top BLAST subjects: 1-aminocyclopropane-1-carboxylate synthase 11 (x10); 1-aminocyclopropane-1-carboxylate synthase 3 (x7); 1-aminocyclopropane-1-carboxylate synthase 4 (x6)
- **function**: Aminotran_1_2 (x10) / 1-aminocyclopropane-1-carboxylate synthase 11 (x10)
- calcification rationale: no direct evidence

### OG0019887  - _LOW_

- sequences: 12; mean length: 297 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: MMR_HSR1(4); Intron_maturas2(4); CysS_C(1)
- top BLAST subjects: Mitochondrial ribosome-associated GTPase 1 (x4); Mitochondrial GTPase 1 (x2)
- **function**: MMR_HSR1 (x4) / Mitochondrial ribosome-associated GTPase 1 (x4)
- calcification rationale: no direct evidence

### OG0019917  - _LOW_

- sequences: 12; mean length: 511 aa
- signal peptide: 8.3%; ≥1 TM region: 25.0% (mean TMRs 0.33)
- top PFAM domains: cNMP_binding(2); Apolipoprotein(1)
- top BLAST subjects: cAMP-dependent protein kinase regulatory subunit (x8)
- **function**: cNMP_binding (x2) / cAMP-dependent protein kinase regulatory subunit (x8)
- calcification rationale: no direct evidence | topology: partial TM (25% TM)

### OG0020622  - _LOW_

- sequences: 11; mean length: 725 aa
- signal peptide: 9.1%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: AAA_30(10); AAA_19(10); UvrD_C_2(10); SH3_13(10); Viral_helicase1(10)
- top BLAST subjects: ATP-dependent RecD2 DNA helicase (x26); RecBCD enzyme subunit RecD (x9); Uncharacterized protein MJ1519 (x2)
- **function**: AAA_30 (x10) / ATP-dependent RecD2 DNA helicase (x26)
- calcification rationale: no direct evidence

### OG0020657  - _LOW_

- sequences: 11; mean length: 339 aa
- signal peptide: 27.3%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence | topology: partially secreted (27% SP)

### OG0020696  - _LOW_

- sequences: 11; mean length: 252 aa
- signal peptide: 9.1%; ≥1 TM region: 27.3% (mean TMRs 0.27)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence | topology: partial TM (27% TM)

### OG0020700  - _LOW_

- sequences: 11; mean length: 387 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: PDZ_2(3); PDZ(2); PDZ_6(2)
- top BLAST subjects: -
- **function**: PDZ_2 (x3)
- calcification rationale: no direct evidence

### OG0020706  - _LOW_

- sequences: 11; mean length: 265 aa
- signal peptide: 0.0%; ≥1 TM region: 45.5% (mean TMRs 2.64)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence | topology: partial TM (45% TM)

### OG0020728  - _LOW_

- sequences: 11; mean length: 299 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0020729  - _LOW_

- sequences: 11; mean length: 282 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0021406  - _LOW_

- sequences: 10; mean length: 230 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: adh_short(10); adh_short_C2(10); KR(9); SDR(5); Epimerase(4)
- top BLAST subjects: 2-dehydro-3-deoxy-D-gluconate 5-dehydrogenase (x17); Gluconate 5-dehydrogenase (x6); 5-keto-D-gluconate 5-reductase (x3)
- **function**: adh_short (x10) / 2-dehydro-3-deoxy-D-gluconate 5-dehydrogenase (x17)
- calcification rationale: no direct evidence

### OG0021517  - _LOW_

- sequences: 10; mean length: 262 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: F-box-like(8); F-box(4)
- top BLAST subjects: F-box only protein 3 (x3); F-box only protein 36 (x2)
- **function**: F-box-like (x8) / F-box only protein 3 (x3)
- calcification rationale: no direct evidence

### OG0021520  - _LOW_

- sequences: 10; mean length: 587 aa
- signal peptide: 0.0%; ≥1 TM region: 10.0% (mean TMRs 1.40)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0021523  - _LOW_

- sequences: 10; mean length: 324 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: PP2C(6)
- top BLAST subjects: Probable protein phosphatase 2C 64 (x5); Probable protein phosphatase 2C 35 (x4); Protein phosphatase 2C and cyclic nucleotide-binding/kinase domain-containing protein (x3)
- **function**: PP2C (x6) / Probable protein phosphatase 2C 64 (x5)
- calcification rationale: no direct evidence

### OG0021531  - _LOW_

- sequences: 10; mean length: 272 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: PDZ_6(10); PDZ_2(7); PDZ(6); PDZ_Tricorn(2); GRASP55_65(2)
- top BLAST subjects: Probable periplasmic serine endoprotease DegP-like (x6); Periplasmic serine endoprotease DegP (x2); Putative serine protease HhoB (x1)
- **function**: PDZ_6 (x10) / Probable periplasmic serine endoprotease DegP-like (x6)
- calcification rationale: no direct evidence

### OG0022355  - _LOW_

- sequences: 9; mean length: 598 aa
- signal peptide: 66.7%; ≥1 TM region: 66.7% (mean TMRs 0.67)
- top PFAM domains: DUF2237(7)
- top BLAST subjects: -
- **function**: DUF2237 (x7)
- calcification rationale: no direct evidence | topology: predominantly secreted (67% SP), membrane-bound (67% TM)

### OG0022381  - _LOW_

- sequences: 9; mean length: 513 aa
- signal peptide: 0.0%; ≥1 TM region: 44.4% (mean TMRs 0.44)
- top PFAM domains: FAD_binding_1(8); NAD_binding_1(8); Flavodoxin_1(7)
- top BLAST subjects: NADPH--cytochrome P450 reductase (x15); NADPH--cytochrome P450 reductase 2 (x8); NADPH--cytochrome P450 reductase 3 (x6)
- **function**: FAD_binding_1 (x8) / NADPH--cytochrome P450 reductase (x15)
- calcification rationale: no direct evidence | topology: partial TM (44% TM)

### OG0022455  - _LOW_

- sequences: 9; mean length: 177 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0022473  - _LOW_

- sequences: 9; mean length: 331 aa
- signal peptide: 33.3%; ≥1 TM region: 33.3% (mean TMRs 0.33)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence | topology: partially secreted (33% SP), partial TM (33% TM)

### OG0022474  - _LOW_

- sequences: 9; mean length: 215 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: F-box-like(8); F-box(6)
- top BLAST subjects: -
- **function**: F-box-like (x8)
- calcification rationale: no direct evidence

### OG0022492  - _LOW_

- sequences: 9; mean length: 472 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Beta-prop_TEP1_2nd(7); WD40(5); Beta-prop_WDR3_1st(5); WD40_Prp19(5); Beta-prop_THOC3(5)
- top BLAST subjects: Myosin heavy chain kinase B (x5); Myosin heavy chain kinase C (x3); F-box/WD repeat-containing protein 7 (x2)
- **function**: Beta-prop_TEP1_2nd (x7) / Myosin heavy chain kinase B (x5)
- calcification rationale: no direct evidence

### OG0022500  - _LOW_

- sequences: 9; mean length: 285 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0022520  - _LOW_

- sequences: 9; mean length: 250 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0022524  - _LOW_

- sequences: 9; mean length: 345 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Pentapeptide(9); Pentapeptide_4(9); Pentapeptide_3(9); Decapeptide(4)
- top BLAST subjects: Uncharacterized protein in mobD 3'region (x71); Uncharacterized protein slr1819 (x33); Uncharacterized protein slr1851 (x23)
- **function**: Pentapeptide (x9) / Uncharacterized protein in mobD 3'region (x71)
- calcification rationale: no direct evidence

### OG0022528  - _LOW_

- sequences: 9; mean length: 320 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0023496  - _LOW_

- sequences: 8; mean length: 411 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: Pkinase(8); PK_Tyr_Ser-Thr(8); Kdo(1)
- top BLAST subjects: Cyclin-G-associated kinase (x12); BMP-2-inducible protein kinase (x3); AP2-associated protein kinase 1 (x3)
- **function**: Pkinase (x8) / Cyclin-G-associated kinase (x12)
- calcification rationale: no direct evidence

### OG0023566  - _LOW_

- sequences: 8; mean length: 475 aa
- signal peptide: 37.5%; ≥1 TM region: 12.5% (mean TMRs 0.12)
- top PFAM domains: RHD3_GTPase(3)
- top BLAST subjects: Protein SEY1 (x5); Protein SEY1 homolog (x4); Protein ROOT HAIR DEFECTIVE 3 homolog 2 (x3)
- **function**: RHD3_GTPase (x3) / Protein SEY1 (x5)
- calcification rationale: no direct evidence | topology: partially secreted (38% SP)

### OG0023594  - _LOW_

- sequences: 8; mean length: 566 aa
- signal peptide: 12.5%; ≥1 TM region: 12.5% (mean TMRs 0.25)
- top PFAM domains: Glyco_transf_8(6); LRR_5(1); Ank_2(1); Ank_4(1); Ank(1)
- top BLAST subjects: Ankyrin repeat, SAM and basic leucine zipper domain-containing protein 1 (x6); Glycosyltransferase 8 domain-containing protein 1 (x4); Probable galacturonosyltransferase-like 4 (x2)
- **function**: Glyco_transf_8 (x6) / Ankyrin repeat, SAM and basic leucine zipper domain-containing protein 1 (x6)
- calcification rationale: no direct evidence

### OG0023657  - _LOW_

- sequences: 8; mean length: 331 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0023790  - _LOW_

- sequences: 8; mean length: 307 aa
- signal peptide: 12.5%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: RCC_reductase(6)
- top BLAST subjects: Red chlorophyll catabolite reductase (Fragment) (x3); Red chlorophyll catabolite reductase 1, chloroplastic (x2)
- **function**: RCC_reductase (x6) / Red chlorophyll catabolite reductase (Fragment) (x3)
- calcification rationale: no direct evidence

### OG0024846  - _LOW_

- sequences: 7; mean length: 333 aa
- signal peptide: 14.3%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: TOD1_MUCI70(4)
- top BLAST subjects: Probable hexosyltransferase MUCI70 (x1)
- **function**: TOD1_MUCI70 (x4) / Probable hexosyltransferase MUCI70 (x1)
- calcification rationale: no direct evidence

### OG0024979  - _LOW_

- sequences: 7; mean length: 492 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0025015  - _LOW_

- sequences: 7; mean length: 167 aa
- signal peptide: 0.0%; ≥1 TM region: 100.0% (mean TMRs 3.14)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence | topology: membrane-bound (100% TM)

### OG0025016  - _LOW_

- sequences: 7; mean length: 217 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: -
- top BLAST subjects: -
- **function**: no detectable homology or domain
- calcification rationale: no direct evidence

### OG0026865  - _LOW_

- sequences: 6; mean length: 492 aa
- signal peptide: 0.0%; ≥1 TM region: 0.0% (mean TMRs 0.00)
- top PFAM domains: JmjC(5); Cupin_8(1)
- top BLAST subjects: Bifunctional arginine demethylase and lysyl-hydroxylase JMJD6 (x7); Bifunctional arginine demethylase and lysyl-hydroxylase PSR (x5); Arginine-specific demethylase JMJ22 (x2)
- **function**: JmjC (x5) / Bifunctional arginine demethylase and lysyl-hydroxylase JMJD6 (x7)
- calcification rationale: no direct evidence

