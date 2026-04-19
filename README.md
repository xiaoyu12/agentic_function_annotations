# Bioinformatics Agent Benchmark for Orthogroup Functional Annotation

This repository contains a controlled comparison of three AI agent workflows for annotating protein orthogroups with complex biological functions and possible relevance to calcification. The benchmark uses the same 73 selected orthogroup FASTA files and the same local bioinformatics evidence for every run, then compares the resulting annotations for biological correctness, missing evidence, hallucinated or over-specific claims, and run-to-run consistency.

The main output is a merged, evidence-audited orthogroup annotation table and a manuscript draft comparing the three agent configurations.

## Benchmark Question

The benchmark evaluates how reliably AI agents can complete the following task:

> The folder `Orthogroups.calcifying_loose_fastas` contains the FASTA files of selected orthogroups. Please annotate the biological functions of each orthogroup and check their relevances to calcification, based on `blast_out`, PFAM domains in `hmm_out`, `deeptmhmm_results` and `signalp_results` in the folder. Generate additional data if necessary.

The three compared methods are:

1. Claude App with Claude Opus 4.7 (`Claude/`)
2. Claude Code CLI with Claude Opus 4.7 and Claude Scientific Skills (`Claude_code/`)
3. Codex App with GPT-5.4 and Claude Scientific Skills (`Codex/`)

Each method has three independent runs.

## Repository Layout

```text
.
|-- README.md
|-- compare_and_merge_annotations.ps1
|-- comparison_curated_final_calls.psv
|-- Orthogroups.calcifying_loose_fastas/
|-- Claude/
|-- Claude_code/
|-- Codex/
|-- comparison_merged_annotation/
`-- biorxiv-paper/
```

Key folders:

- `Orthogroups.calcifying_loose_fastas/`: shared FASTA and evidence inputs used by all agent runs.
- `Claude/`: outputs from Claude App Opus 4.7 runs.
- `Claude_code/`: outputs from Claude Code CLI Opus 4.7 runs with Claude Scientific Skills.
- `Codex/`: outputs from Codex App GPT-5.4 runs with Claude Scientific Skills.
- `comparison_merged_annotation/`: normalized comparison tables, recomputed evidence summaries, and final merged annotations.
- `biorxiv-paper/`: LaTeX manuscript, bibliography, generated figures, supplement tables, and manuscript PDF.

## Main Outputs

The most important benchmark outputs are in `comparison_merged_annotation/`:

- `final_orthogroup_annotations.tsv`: final curated annotation table for all 73 orthogroups.
- `model_run_comparison.tsv`: per-orthogroup comparison across all nine agent runs.
- `method_consistency_summary.tsv`: run-level label distributions and score summaries.
- `recomputed_evidence.tsv`: evidence re-aggregation from FASTA, BLAST, Pfam/HMMER, SignalP, and DeepTMHMM.
- `comparison_and_final_annotation_report.md`: human-readable audit report.

The manuscript is in `biorxiv-paper/`:

- `main.tex`: manuscript source.
- `references.bib`: bibliography.
- `main.pdf`: compiled manuscript.
- `figures/`: generated manuscript figures.
- `supplement/`: generated supplementary LaTeX tables.
- `sources/web_research_20260418.md`: web and literature sources used while drafting the paper.

## Summary of Findings

All nine agent runs covered all 73 orthogroups, so the main differences were interpretive rather than file-retrieval failures. The final merged annotation classified:

- High relevance: 3 orthogroups
- Moderate relevance: 9 orthogroups
- Watchlist: 18 orthogroups
- Low relevance: 43 orthogroups

The strongest high-confidence candidates were sulfatase or sulfotransferase families and an FG-GAP/integrin-like surface candidate. Common overclaim patterns included:

- Treating pentapeptide-repeat proteins as direct calcification candidates without supporting evidence.
- Promoting secreted housekeeping enzymes based on SignalP alone.
- Over-interpreting single weak collagen or adhesive BLAST hits.
- Transferring low-complexity or giant-protein BLAST labels without domain support.

The skill-enabled coding-agent workflows improved evidence handling, reproducible scripts, and auditability, but they did not eliminate biological overinterpretation.

## Reproducing the Comparison

Run the comparison and merge workflow from the repository root:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\compare_and_merge_annotations.ps1
```

This regenerates the comparison outputs in `comparison_merged_annotation/`.

## Regenerating Manuscript Figures and Supplementary Tables

The manuscript assets are generated from the comparison TSV files using PowerShell and .NET drawing utilities:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\biorxiv-paper\scripts\generate_manuscript_assets.ps1
```

This regenerates:

- `biorxiv-paper/figures/final_relevance_distribution.png`
- `biorxiv-paper/figures/run_relevance_distribution.png`
- `biorxiv-paper/figures/within_method_agreement.png`
- `biorxiv-paper/figures/candidate_evidence_heatmap.png`
- `biorxiv-paper/supplement/final_annotations_table.tex`
- `biorxiv-paper/supplement/method_summary_table.tex`
- `biorxiv-paper/supplement/pairwise_agreement_table.tex`

## Building the Manuscript

From `biorxiv-paper/`, compile the LaTeX manuscript with:

```powershell
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The compiled PDF is written to:

```text
biorxiv-paper/main.pdf
```

## Requirements

The comparison and asset-generation scripts are designed for Windows PowerShell. Building the manuscript requires a LaTeX distribution such as MiKTeX with `pdflatex` and `bibtex` available.

No Python or R environment is required for the current figure-generation workflow.

## Citation

If you use this repository, cite the manuscript in `biorxiv-paper/main.pdf` or the repository itself until a preprint DOI is available.

## Notes for Reuse

This benchmark is intended as a reproducible audit of AI-assisted biological annotation, not as a claim that any single agent configuration is generally superior. For similar protein annotation tasks, use deterministic evidence tables, explicit relevance-tier definitions, multiple independent agent runs, and expert review of high-impact biological claims.
