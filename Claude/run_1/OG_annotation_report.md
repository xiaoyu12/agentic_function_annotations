# Orthogroup functional annotation — calcifying-loose set

## Context

Seventy-three orthogroups from the "calcifying_loose" selection were annotated using:
- **BLAST** against UniProt/Swiss-Prot (`blast_out/*.blast.tsv`)
- **PFAM** domain HMM search (`hmm_out/*.pfam.tbl`)
- **SignalP-6.0** signal-peptide predictions (`signalp_results/*/output.gff3`)
- **DeepTMHMM** transmembrane topology (`deeptmhmm_results/*.TMRs.gff3`)

Evidence was aggregated per OG with `aggregate_annotations.py` into `annotation_summary.tsv` / `.json` (included in this folder).

### Taxon coding (22 haptophytes in total)

| Calcifying (coccolithophores) | Non-calcifying haptophytes |
|---|---|
| Callep1130_1 *(Calcidiscus leptoporus)* | Chrlea393_1, Chrsim136_1 *(Chrysochromulina spp.)* |
| Emihu1 *(Emiliania huxleyi)* | Imarot2298_1, Imarot4477_3, Imarot4780_1, Imarot704_1 *(Imantonia rotunda)* |
| Gepoce1 *(Gephyrocapsa oceanica)* | IsochDm2_1, Isogal1 *(Isochrysis spp.)* |
| Hymcor1016_1 *(Hymenomonas coronata)* | Pavlova2550_1, Pavlova5622_1, Pavlova6257_1 *(Pavlovales)* |
| Ochro3194_1 *(Ochrosphaera neapolitana)* | Phaant1, Phaglo1 *(Phaeocystis spp.)* |
| Pleelo874_1 *(Pleurochrysis elongata)* | Platy1217_1 *(Platychrysis)* |
| Umbfol2878_1 *(Umbellosphaera)* | Psesor5268_1 *(Pseudoisochrysis)* |

### Calcification-relevance score legend

- **HIGH** — function plausibly and directly involved in coccolith biomineralisation (Ca²⁺/HCO₃⁻ transport, matrix acidic polysaccharide synthesis, coccolith-vesicle / exocytosis machinery, known coccolith-associated protein families).
- **MEDIUM** — function could support calcification indirectly (redox/pH homeostasis, secreted glycan-modifying enzymes, cell-wall / extracellular matrix remodelling, vesicle trafficking) but is not calcification-specific.
- **LOW** — common housekeeping, nucleic-acid, general metabolism or cytoskeletal roles with no obvious link to biomineralisation.
- **UNKNOWN** — no reliable functional signal.
- **ARTEFACT?** — profile consistent with contamination / repeat element / prokaryotic transfer.

---

## Per-orthogroup annotations

### OG0000049 — 502 seq · 7 spp (Pleelo 313, Callep 124, Umbfol 42, Gepoce 10, Ochro 8, Hymcor 4, Emihu 1)
- **BLAST cov 4 %**, scattered best-hits (Endonuclease III, FKBP15, Reticulocyte-binding protein); **PFAM cov 3 %** (FKBP15, CCDC39, ODF2 coiled-coil, MAD).
- **SignalP 0 %**, **TM 1 %** (mostly intracellular), no domain dominates.
- **Interpretation:** heterogeneous cluster enriched in coiled-coil / rapidly evolving sequences, highly amplified in *Pleurochrysis* and *Calcidiscus* — likely lineage-specific / repetitive proteins. Relevance: **UNKNOWN** (possibly ARTEFACT? from rapidly expanded low-complexity family).

### OG0001332 — 91 seq · 8 spp (Callep 76 dominant)
- Best hit **Warthog-4 / WRT1 / Hedgehog-like**, PFAM **Hint domain** (PF01079, 85/91 seqs).
- **84 % TM**, 1 signal peptide.
- **Interpretation:** Hint-domain / Hedgehog-like self-splicing transmembrane proteins, massively expanded in *Calcidiscus*. Could participate in cell-surface patterning but no direct calcification link.
- Relevance: **LOW–MEDIUM** (surface-anchored, haptophyte-specific expansion).

### OG0001976 — 73 seq · 7 spp (all calcifying)
- Pentapeptide-repeat proteins (PF00805, PF13599, PF13576) and Decapeptide repeats; best hits to cyanobacterial *Synechocystis* Y1819/Y1152 and *Pantoea* YMO3.
- **SignalP 5 %**, **TM 4 %**, cytosolic.
- **Interpretation:** Pentapeptide-repeat proteins (β-helical solenoids) — in cyanobacteria associated with calcium and carbonate binding surfaces; in haptophytes reported in EST libraries of calcifying cells.
- Relevance: **MEDIUM–HIGH** (same family is a well-known candidate for Ca²⁺/carbonate interaction in cyanobacterial calcification; expansion here is restricted to coccolithophores).

### OG0002887 — 59 seq · 9 spp
- Best hit **PIF1 DNA helicase**, domain **PIF1** + **Helitron** helicase.
- 0 % TM, 2 % SignalP.
- **Interpretation:** nuclear DNA helicase / Helitron transposase-like.
- Relevance: **LOW**.

### OG0003717 — 51 seq · 8 spp
- No BLAST or PFAM hits; 0 SignalP; 8 % TM.
- **Interpretation:** haptophyte-specific orphans with weak membrane signal. Cannot assign.
- Relevance: **UNKNOWN** (candidate taxon-restricted gene worth follow-up).

### OG0009246 — 29 seq · 7 spp
- Best BLAST **Adhesive plaque matrix protein FP1 (mussel foot)**; sparse PFAM (tRNA-synt_2 noise).
- **31 % SignalP**, **14 % TM** — suggests secreted/acidic protein.
- **Interpretation:** secreted Gly/Ala/Ser-rich adhesive-plaque-like proteins; repeat-rich signature characteristic of biomineralisation matrix proteins.
- Relevance: **HIGH** candidate coccolith-matrix / extracellular adhesive glycoprotein.

### OG0009301 — 29 seq · 7 spp
- Minimal BLAST (spurious RPB1), PFAM **SHOCT (Short C-terminal)** + **PRIMA1 (proline-rich membrane anchor)**.
- **72 % SignalP**, **52 % TM** — a highly secreted / membrane-anchored family.
- **Interpretation:** small secreted proline/SHOCT-anchor proteins, often associated with cell-surface peripheral membrane complexes.
- Relevance: **MEDIUM–HIGH** (secreted membrane-anchored protein family restricted to calcifying spp. fits coccolith-vesicle / cell-surface calcification scaffold profile).

### OG0009816 — 28 seq · 8 spp
- Best hit **MNN14 α-1,3-mannosyltransferase**, PFAM **Mannosyl_trans3 (PF11051, 15/28)**.
- **7 % SignalP**, **21 % TM** — Golgi-anchored glycosyltransferase.
- **Interpretation:** α-1,3-mannosyltransferase.
- Relevance: **MEDIUM** (N/O-glycosylation of secreted coccolith proteins is a plausible calcification support role).

### OG0010177 — 28 seq · 7 spp (21/28 in Chrlea!)
- Dense **ankyrin-repeat** protein (PF00023/PF12796/…); best hit *Rickettsia* ankyrin repeat / ANR17_HUMAN.
- 0 % SignalP, 0 % TM.
- **Interpretation:** cytosolic ankyrin-repeat scaffolds; dominant in *Chrysochromulina* suggests this OG is mostly non-calcifier-driven.
- Relevance: **LOW**.

### OG0010867 — 27 seq · 7 spp (all calcifying)
- **Pentapeptide / Decapeptide repeat proteins**, 100 % PFAM cov, 100 % BLAST cov (Y1819_SYNY3, YMO3_PANSE).
- 0 TM, 4 % SignalP.
- **Interpretation:** second pentapeptide-repeat family (companion to OG0001976). Restricted to coccolithophores.
- Relevance: **MEDIUM–HIGH** (see OG0001976).

### OG0010991 — 27 seq · 7 spp
- No BLAST / PFAM, 15 % SignalP, 4 % TM.
- **Interpretation:** unannotated, weak secretion signal.
- Relevance: **UNKNOWN** (taxon-restricted candidate).

### OG0011061 — 27 seq · 9 spp
- Best hit PWWP4 (human) + **Exostosin GT47** glycosyltransferase PFAM.
- **33 % SignalP**, **100 % TM** (every sequence has a TM span).
- **Interpretation:** GT47 **xylosyl/glucuronyltransferase** integral-membrane enzyme (Golgi type II) — builds heparan-sulfate-like / pectic side-chain polysaccharides.
- Relevance: **HIGH** (coccolith-associated polysaccharides [CAP] are sulphated / uronic-acid glycans; GT47 is a strong biosynthesis candidate).

### OG0011197 — 27 seq · 8 spp (Hymcor 14 dominant)
- Best hit **RAP/2-OG-Fe(II) oxygenase** (PF03171).
- 4 % SignalP, 7 % TM.
- **Interpretation:** Fe(II)/2-oxoglutarate-dependent dioxygenase, often chloroplastic.
- Relevance: **LOW** (metabolic, unlikely calcification).

### OG0013037 — 24 seq · 6 spp
- Weak titin-like coiled-coil hits; PFAM **Hydrolase_4 (α/β-serine-aminopeptidase)**.
- 21 % SignalP, 25 % TM.
- **Interpretation:** serine peptidase with partial secretion.
- Relevance: **LOW–MEDIUM**.

### OG0014155 — 22 seq · 8 spp
- Best hit **Nucleoside hydrolase NSH3/4/5**; PFAM **IU_nuc_hydro (PF01156, 17/22)**.
- **50 % SignalP**, 9 % TM.
- **Interpretation:** secreted inosine-uridine nucleoside hydrolase (purine salvage).
- Relevance: **LOW**.

### OG0014250 — 22 seq · 8 spp
- Best hit **Glutamyl-tRNA(Gln) amidotransferase A (GatA)**; PFAM **Amidase (PF01425, 17/22)**.
- **68 % SignalP**, 5 % TM.
- **Interpretation:** secreted amidase (nitrilase-superfamily); GatA hit is likely spurious — catalytic PFAM "Amidase" is the cue. Extracellular organic-N / amide hydrolysis.
- Relevance: **LOW–MEDIUM** (secreted but not clearly calcification-linked).

### OG0015153 — 20 seq · 9 spp
- Poor BLAST (Piccolo/HAC6 low-complexity); no PFAM.
- **85 % TM**, 5 % SignalP.
- **Interpretation:** multi-pass membrane protein of unknown function, broadly distributed.
- Relevance: **UNKNOWN** (integral-membrane candidate worth investigation).

### OG0016203 — 18 seq · 9 spp
- Best hit **Collagen α-2(IV)**; PFAM **Collagen triple-helix** + **RCC1** repeats.
- 6 % SignalP, **83 % TM**.
- **Interpretation:** collagen-like glycine-rich repeat proteins, membrane-anchored; collagen repeats are structural scaffolds.
- Relevance: **MEDIUM–HIGH** (extracellular structural collagen-like scaffolds are documented in some biomineralising protists; intriguing candidate).

### OG0016211 — 18 seq · 9 spp
- Minimal PFAM (**SUR7/PalI family**); 11 % SignalP, **100 % TM**.
- **Interpretation:** SUR7/PalI multi-pass membrane microdomain proteins (eisosome / MCC-equivalent in yeast) — lipid-raft scaffolds in plasma membrane.
- Relevance: **MEDIUM** (plasma-membrane domain organisation could be relevant to coccolith-vesicle docking but speculative).

### OG0017138 — 16 seq · 8 spp
- Best hit **Arylsulfatase B (ARSB)**, PFAM **Sulfatase (PF00884, 14/16)**.
- 12 % SignalP, 0 TM.
- **Interpretation:** arylsulfatase (hydrolyses sulfate esters).
- Relevance: **MEDIUM** (sulfated coccolith-associated polysaccharides exist; sulfatases regulate turnover/maturation of sulfated glycans).

### OG0017238 — 16 seq · 9 spp
- Best hit **tRNA ligase 1 (RNL_ARATH)**, no PFAM.
- 0 % SignalP / TM.
- **Interpretation:** intracellular RNA ligase.
- Relevance: **LOW**.

### OG0017305 — 16 seq · 8 spp
- PFAM **Stf0 sulfotransferase (PF09037, 3/16)**; no BLAST.
- 6 % SignalP, 0 TM.
- **Interpretation:** sulfotransferase (adds sulfate to sugars/metabolites).
- Relevance: **MEDIUM–HIGH** (sulfation of coccolith-associated polysaccharides is a known calcification-linked modification).

### OG0017361 — 16 seq · 8 spp
- Best hit **CdhR / InhA transcriptional regulator**; PFAM **DJ-1/PfpI (PF01965, 16/16)** + HTH_AraC.
- 6 % SignalP, 12 % TM.
- **Interpretation:** DJ-1/PfpI-family cysteine peptidase / glyoxalase (redox-sensing, protein deglycase). Broad detoxification / stress role.
- Relevance: **LOW–MEDIUM**.

### OG0017362 — 16 seq · 9 spp
- Best hit APE1L endonuclease + **RCC_reductase (red chlorophyll catabolite reductase)**.
- **44 % SignalP**, 12 % TM.
- **Interpretation:** mixed annotation — likely chlorophyll catabolism / senescence pathway, possibly also a secreted subfamily.
- Relevance: **LOW**.

### OG0017884 — 15 seq · 7 spp
- PFAM **GCN5-related N-acetyltransferase (Acetyltransf_1/7)**; best hit mycothiol acetyltransferase.
- 0 SignalP, 27 % TM.
- **Interpretation:** GNAT acetyltransferase.
- Relevance: **LOW**.

### OG0017903 — 15 seq · 8 spp
- **Flavohemoprotein (HMP)**, PFAM **FAD_binding_6 + Globin + NAD_binding_1**.
- 0 SignalP / TM.
- **Interpretation:** bacterial-type flavohemoglobin — NO detoxification / redox.
- Relevance: **LOW**.

### OG0017965 — 15 seq · 6 spp
- PFAM **Kelch / NANM** β-propeller, best hit KEAP1 / NANM (N-acetylneuraminate epimerase).
- 0 SignalP / TM.
- **Interpretation:** Kelch-β-propeller redox/substrate-adaptor proteins.
- Relevance: **LOW**.

### OG0018519 — 14 seq · 6 spp
- No BLAST / PFAM; **93 % TM**.
- **Interpretation:** uncharacterised multi-pass membrane protein; coccolith-lineage-only.
- Relevance: **UNKNOWN** (strong membrane candidate — flagged for experimental follow-up).

### OG0018986 — 13 seq · 7 spp
- Best hit **metabotropic glutamate receptor-like E (GRLE)**; PFAM **FG-GAP + FG-GAP_3 + 7tm_3** (GPCR) + Insecticide toxin mid.
- **54 % SignalP**, **77 % TM**.
- **Interpretation:** 7TM GPCR-like receptor with FG-GAP (integrin) repeats — extracellular-sensing receptor family.
- Relevance: **MEDIUM** (calcium/pH sensing at plasma membrane plausible).

### OG0019067 — 13 seq · 8 spp
- **LRR** proteins (PF13516), best hit **NLRC3** (NLR innate-immunity LRR).
- 0 SignalP, 38 % TM.
- **Interpretation:** intracellular leucine-rich-repeat protein.
- Relevance: **LOW**.

### OG0019078 — 13 seq · 9 spp
- Best hit **Elongation factor 1-γ**, PFAM **GST_N + GST_C**.
- 0 SignalP / TM.
- **Interpretation:** EF1γ (GST-fold EF1B subunit). Housekeeping.
- Relevance: **LOW**.

### OG0019174 — 13 seq · 7 spp
- Best hit **GatA (Glu-tRNA(Gln) amidotransferase A)**; PFAM **Amidase (PF01425, 13/13)**.
- **38 % SignalP**, 0 TM.
- **Interpretation:** amidase / nitrilase-superfamily; partial secretion.
- Relevance: **LOW**.

### OG0019180 — 13 seq · 7 spp
- No BLAST / PFAM; no SignalP / TM.
- **Interpretation:** orphan cytosolic family.
- Relevance: **UNKNOWN**.

### OG0019790 — 12 seq · 9 spp
- Best hit **HECD1 HECT-type E3 ubiquitin ligase**; PFAM **HECT + Ankyrin repeats**.
- 0 SignalP / TM.
- **Interpretation:** ubiquitin-protein ligase.
- Relevance: **LOW**.

### OG0019817 — 12 seq · 8 spp
- Best hit **Aurachin-C monooxygenase (AuaG)**; PFAM **FAD_binding_3 + NAD_binding_8**.
- **67 % SignalP**, 8 % TM.
- **Interpretation:** FAD/NAD(P)-dependent aromatic monooxygenase / hydroxylase, secreted.
- Relevance: **LOW–MEDIUM**.

### OG0019873 — 12 seq · 9 spp
- Best hit **1-aminocyclopropane-1-carboxylate synthase (ACS)**; PFAM **Aminotran_1_2**.
- 0 SignalP / TM.
- **Interpretation:** PLP-dependent aminotransferase / ACC-synthase homolog.
- Relevance: **LOW**.

### OG0019887 — 12 seq · 8 spp
- Best hit **Mitochondrial GTPase 1 (MTG1)**; PFAM **MMR_HSR1** + **Intron_maturas2**.
- 0 SignalP / TM.
- **Interpretation:** ribosome-biogenesis GTPase / group-II intron maturase (organellar).
- Relevance: **LOW**.

### OG0019917 — 12 seq · 6 spp
- Best hit **cAMP-dependent PK regulatory subunit (KAPR)**; PFAM **cNMP_binding + Apolipoprotein**.
- 8 % SignalP, 25 % TM.
- **Interpretation:** PKA regulatory subunit (cyclic-nucleotide binding).
- Relevance: **LOW**.

### OG0020622 — 11 seq · 8 spp
- Best hit **RecD2 ATP-dependent DNA helicase**, PFAM **AAA_30 + UvrD_C_2 + Viral_helicase1**.
- 9 % SignalP, 0 TM.
- **Interpretation:** SF1 DNA helicase.
- Relevance: **LOW**.

### OG0020657 — 11 seq · 8 spp
- No BLAST / PFAM. 27 % SignalP, 0 TM.
- **Interpretation:** partially secreted orphan.
- Relevance: **UNKNOWN**.

### OG0020696 — 11 seq · 8 spp
- No BLAST / PFAM. 9 % SignalP, 27 % TM.
- **Interpretation:** membrane orphan family.
- Relevance: **UNKNOWN**.

### OG0020700 — 11 seq · 8 spp
- No BLAST; PFAM **PDZ (×3)**.
- 0 SignalP / TM.
- **Interpretation:** cytosolic PDZ-scaffold.
- Relevance: **LOW**.

### OG0020703 — 11 seq · 8 spp
- Best hit **Heparan-sulfate 6-O-sulfotransferase**, PFAM **Sulfotransfer_1/3**.
- **36 % SignalP**, 0 TM.
- **Interpretation:** Golgi-lumen heparan-sulfate 6-O-sulfotransferase.
- Relevance: **MEDIUM–HIGH** (sulfation of CAP-type acidic polysaccharides).

### OG0020706 — 11 seq · 9 spp
- No BLAST / PFAM; 0 SignalP, **45 % TM**.
- **Interpretation:** multi-pass membrane orphan.
- Relevance: **UNKNOWN**.

### OG0020728 — 11 seq · 8 spp
- No BLAST / PFAM / SignalP / TM.
- **Interpretation:** cytosolic orphan.
- Relevance: **UNKNOWN**.

### OG0020729 — 11 seq · 7 spp
- No BLAST / PFAM / SignalP / TM.
- **Interpretation:** cytosolic orphan.
- Relevance: **UNKNOWN**.

### OG0021347 — 10 seq · 8 spp
- Best hit **Trypsin 3A1**; PFAM **Trypsin (PF00089, 10/10)**.
- **50 % SignalP**, 0 TM.
- **Interpretation:** secreted serine (trypsin) protease.
- Relevance: **LOW–MEDIUM** (extracellular proteolysis can modulate biomineral matrix).

### OG0021406 — 10 seq · 8 spp
- Best hit **2-dehydro-3-deoxy-D-gluconate 5-dehydrogenase (KduD)**; PFAM **adh_short + KR**.
- 0 SignalP / TM.
- **Interpretation:** SDR superfamily dehydrogenase (sugar-acid metabolism).
- Relevance: **LOW**.

### OG0021517 — 10 seq · 8 spp
- Best hit **FBX3**; PFAM **F-box** (8/10).
- 0 SignalP / TM.
- **Interpretation:** F-box protein (ubiquitin ligase substrate adaptor).
- Relevance: **LOW**.

### OG0021520 — 10 seq · 8 spp
- No BLAST / PFAM; 0 SignalP, 10 % TM.
- **Interpretation:** orphan.
- Relevance: **UNKNOWN**.

### OG0021523 — 10 seq · 7 spp
- Best hit **PP2C protein phosphatase**; PFAM **PP2C (PF00481, 6/10)**.
- 0 SignalP / TM.
- **Interpretation:** Mg²⁺/Mn²⁺-dependent Ser/Thr phosphatase.
- Relevance: **LOW**.

### OG0021531 — 10 seq · 6 spp
- Best hit **ZO-3 tight-junction** (PDZ); PFAM **PDZ_6/PDZ_2/PDZ/PDZ_Tricorn/GRASP55_65 (10/10)**.
- 0 SignalP / TM.
- **Interpretation:** multi-PDZ scaffold (tricorn-protease / GRASP-like).
- Relevance: **LOW–MEDIUM** (PDZ scaffolds organise membrane-channel complexes; speculative).

### OG0022355 — 9 seq · 8 spp
- No BLAST; PFAM **DUF2237** (7/9).
- **67 % SignalP**, 67 % TM.
- **Interpretation:** DUF2237 is a family of secreted/membrane proteins of unknown function (often γ-glutamyl-cyclotransferase-like).
- Relevance: **MEDIUM** (strong secretion + TM + bacterial-like fold).

### OG0022381 — 9 seq · 8 spp
- Best hit **NADPH-cytochrome P450 reductase**; PFAM **FAD_binding_1 + NAD_binding_1 + Flavodoxin_1**.
- 0 SignalP, 44 % TM.
- **Interpretation:** ER-membrane cytochrome-P450 reductase.
- Relevance: **LOW**.

### OG0022455 — 9 seq · 8 spp
- No evidence — complete orphan.
- Relevance: **UNKNOWN**.

### OG0022473 — 9 seq · 7 spp
- No BLAST / PFAM; 33 % SignalP, 33 % TM.
- **Interpretation:** uncharacterised secreted/membrane family.
- Relevance: **UNKNOWN** (worth further examination).

### OG0022474 — 9 seq · 8 spp
- No BLAST; PFAM **F-box / F-box-like**.
- 0 SignalP / TM.
- **Interpretation:** F-box adaptor.
- Relevance: **LOW**.

### OG0022492 — 9 seq · 8 spp
- Best hit **MHCKB (Myosin heavy-chain kinase B)** — α-kinase/WD40 scaffold; PFAM **WD40 / β-propeller** ×5.
- 0 SignalP / TM.
- **Interpretation:** WD40-repeat α-kinase scaffold.
- Relevance: **LOW**.

### OG0022500 — 9 seq · 9 spp
- No evidence — orphan 1-per-species.
- Relevance: **UNKNOWN** (single-copy, pan-haptophyte orphan — possibly conserved coccolith-lineage protein; flagged).

### OG0022520 — 9 seq · 9 spp
- No evidence — orphan 1-per-species.
- Relevance: **UNKNOWN** (see OG0022500).

### OG0022524 — 9 seq · 7 spp
- Pentapeptide-repeat family (PF00805, 9/9). Best hit YMO3.
- 0 SignalP / TM.
- **Interpretation:** third pentapeptide-repeat paralog subfamily.
- Relevance: **MEDIUM–HIGH** (see OG0001976 / OG0010867).

### OG0022528 — 9 seq · 7 spp
- No BLAST / PFAM / SignalP / TM.
- **Interpretation:** orphan.
- Relevance: **UNKNOWN**.

### OG0023496 — 8 seq · 7 spp
- Best hit **BMP2K / GAK cyclin-G-associated kinase**; PFAM **Pkinase (8/8)**.
- 0 SignalP / TM.
- **Interpretation:** Ser/Thr protein kinase (GAK/BMP2K subfamily — clathrin-uncoating kinase).
- Relevance: **LOW–MEDIUM** (clathrin-uncoating kinases participate in vesicle trafficking; relevant if coccolith-vesicle dynamics require clathrin-mediated steps, but speculative).

### OG0023566 — 8 seq · 8 spp
- Best hit **SEY1 atlastin-like GTPase**; PFAM **RHD3_GTPase**.
- 38 % SignalP, 12 % TM.
- **Interpretation:** Atlastin / Sey1-type ER-membrane fusion GTPase.
- Relevance: **LOW**.

### OG0023594 — 8 seq · 6 spp
- Best hit **ASZ1 (Ankyrin/SAM)**; PFAM **Glyco_transf_8 (PF01501, 6/8)** + Ankyrin.
- 12 % SignalP, 12 % TM.
- **Interpretation:** GT8 glycosyltransferase (galactosyl/glucosyl-transferase family; includes LgtC, GolS, pectin-rhamnogalacturonan transferases).
- Relevance: **MEDIUM** (GT8 family includes enzymes that make uronic-acid-rich polysaccharides — relevant to CAP biosynthesis).

### OG0023657 — 8 seq · 8 spp
- No evidence.
- Relevance: **UNKNOWN**.

### OG0023790 — 8 seq · 8 spp
- Best hit **Red chlorophyll catabolite reductase (RCCR)**; PFAM **RCC_reductase (6/8)**.
- 12 % SignalP, 0 TM.
- **Interpretation:** chlorophyll-catabolism reductase (plastidial).
- Relevance: **LOW**.

### OG0024846 — 7 seq · 7 spp
- Best hit **MUC70 (Arabidopsis MUCI70 hexosyltransferase)**; PFAM **TOD1_MUCI70 (PF04765)**.
- 14 % SignalP, 0 TM.
- **Interpretation:** TOD1/MUCI70-family glycosyltransferase (pectic/mucilage polysaccharide biosynthesis).
- Relevance: **MEDIUM** (candidate CAP-polysaccharide biosynthesis enzyme; the Arabidopsis homolog is a seed-mucilage hexosyltransferase).

### OG0024979 — 7 seq · 7 spp
- No evidence.
- Relevance: **UNKNOWN**.

### OG0025015 — 7 seq · 7 spp
- No BLAST / PFAM; **100 % TM**; 0 SignalP.
- **Interpretation:** single-copy-per-species multi-pass membrane protein — 7 species, 1 seq each.
- Relevance: **UNKNOWN — HIGH PRIORITY** (strict single-copy membrane orphan shared by all sampled calcifiers is the exact topology expected for a coccolith-vesicle integral transporter; flag for experimental follow-up).

### OG0025016 — 7 seq · 7 spp
- No evidence.
- Relevance: **UNKNOWN** (single-copy orphan; lower priority than OG0025015 because no TM).

### OG0025059 — 7 seq · 6 spp
- No BLAST / PFAM; **71 % SignalP** (5/7 have high-confidence signal peptide); 0 TM.
- **Interpretation:** small secreted orphan protein shared across calcifiers.
- Relevance: **HIGH candidate** (secreted, calcifier-restricted, taxon-specific — matches the profile of a novel coccolith-matrix protein; highly recommended for experimental follow-up).

### OG0026865 — 6 seq · 6 spp
- Best hit **JmjD6 arginine demethylase / lysyl-hydroxylase**; PFAM **JmjC + Cupin_8**.
- 0 SignalP / TM.
- **Interpretation:** JmjC-domain 2-OG-Fe(II) hydroxylase/demethylase.
- Relevance: **LOW**.

---

## Calcification-priority shortlist

| Priority | OG | Why |
|---|---|---|
| ★★★ | **OG0025059** | Small **secreted (71 %) calcifier-restricted orphan** — classic novel coccolith-matrix candidate. |
| ★★★ | **OG0025015** | Single-copy, 100 %-TM orphan conserved across 7 calcifiers — putative coccolith-vesicle transporter. |
| ★★★ | **OG0011061** | GT47 + 100 % TM + 33 % SignalP — Golgi glycosyltransferase for sulfated/uronic coccolith polysaccharides. |
| ★★ | **OG0009246** | Repeat-rich, 31 % SignalP adhesive-plaque-like family — coccolith-matrix glycoprotein profile. |
| ★★ | **OG0009301** | 72 % SignalP, 52 % TM secreted/anchored family (SHOCT/PRIMA1). |
| ★★ | **OG0001976 / OG0010867 / OG0022524** | Three pentapeptide-repeat paralog groups, **calcifier-specific**; Ca²⁺/carbonate-surface candidates. |
| ★★ | **OG0017305 / OG0020703 / OG0017138** | Sulfotransferases / arylsulfatase — CAP sulfation enzymes. |
| ★★ | **OG0023594 / OG0024846 / OG0009816** | GT8 / GT-MUCI70 / GT-mannosyl — acidic-polysaccharide biosynthesis candidates. |
| ★★ | **OG0016203** | Collagen-like Gly-rich + RCC1 repeats, 83 % TM — potential extracellular scaffold. |
| ★ | **OG0018986** | 7TM-GPCR + FG-GAP integrin repeats — extracellular Ca²⁺/pH sensor. |
| ★ | **OG0018519 / OG0015153 / OG0016211 / OG0020706 / OG0022473** | Uncharacterised membrane orphans restricted to coccolithophores — worth targeted follow-up. |

---

## Outputs in this folder

- `aggregate_annotations.py` — script that parses BLAST/PFAM/SignalP/DeepTMHMM outputs and produces summaries.
- `annotation_summary.tsv` — one row per OG with top BLAST hit, top 5 PFAM domains, SignalP and TM counts, and species composition.
- `annotation_summary.json` — same data as a machine-readable JSON list.
- `OG_annotation_report.md` — this narrative report.
