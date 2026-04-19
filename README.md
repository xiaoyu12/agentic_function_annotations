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
|-- Orthogroups.calcifying_loose_fastas/
|-- Claude/
|-- Claude_code/
|-- Codex/
|-- comparison_merged_annotation/
```

Key folders:

- `Orthogroups.calcifying_loose_fastas/`: shared FASTA and evidence inputs used by all agent runs.
- `Claude/`: outputs from Claude App Opus 4.7 runs.
- `Claude_code/`: outputs from Claude Code CLI Opus 4.7 runs with Claude Scientific Skills.
- `Codex/`: outputs from Codex App GPT-5.4 runs with Claude Scientific Skills.
- `comparison_merged_annotation/`: normalized comparison tables, recomputed evidence summaries, and final merged annotations.


## Main Outputs

The most important benchmark outputs are in `comparison_merged_annotation/`:

- `final_orthogroup_annotations.tsv`: final curated annotation table for all 73 orthogroups.
- `model_run_comparison.tsv`: per-orthogroup comparison across all nine agent runs.
- `method_consistency_summary.tsv`: run-level label distributions and score summaries.
- `recomputed_evidence.tsv`: evidence re-aggregation from FASTA, BLAST, Pfam/HMMER, SignalP, and DeepTMHMM.
- `comparison_and_final_annotation_report.md`: human-readable audit report.
- `dedman_review/dedman_paper_annotation_review.md`: literature-based update after reviewing Dedman et al. (2024) coccolith matrix proteomics.

## Summary of Findings

All nine agent runs covered all 73 orthogroups, so the main differences were interpretive rather than file-retrieval failures. The final merged annotation classified:

- High relevance: 3 orthogroups
- Moderate relevance: 9 orthogroups
- Watchlist: 21 orthogroups
- Low relevance: 40 orthogroups

The strongest high-confidence candidates were sulfatase or sulfotransferase families and an FG-GAP/integrin-like surface candidate. Common overclaim patterns included:

- Treating pentapeptide-repeat proteins as direct calcification candidates. After reviewing Dedman et al. (2024), these are best handled as watchlist candidates because the motif has coccolith-matrix support, but orthogroup-specific function remains unproven.
- Promoting secreted housekeeping enzymes based on SignalP alone.
- Over-interpreting single weak collagen or adhesive BLAST hits.
- Transferring low-complexity or giant-protein BLAST labels without domain support.

The skill-enabled coding-agent workflows improved evidence handling, reproducible scripts, and auditability, but they did not eliminate biological overinterpretation.


## Reproducing the Merge

The deterministic merge workflow is available in both PowerShell and Python.

PowerShell:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\compare_and_merge_annotations.ps1
```

Python 3.13+:

```powershell
& 'C:\Users\zhang\AppData\Local\Programs\Python\Python313\python.exe' .\merge_best_evidence.py
```

To write a trial run to a separate directory:

```powershell
& 'C:\Users\zhang\AppData\Local\Programs\Python\Python313\python.exe' .\merge_best_evidence.py --out-dir .\comparison_merged_annotation_python_check
```

The Python script reads `comparison_curated_final_calls.psv`, recomputes local FASTA/BLAST/Pfam/SignalP/DeepTMHMM evidence, parses the nine agent outputs, applies the Dedman et al. (2024) literature update layer, and regenerates the final annotation table and audit report.


## Notes for Reuse

This benchmark is intended as a reproducible audit of AI-assisted biological annotation, not as a claim that any single agent configuration is generally superior. For similar protein annotation tasks, use deterministic evidence tables, explicit relevance-tier definitions, multiple independent agent runs, and expert review of high-impact biological claims.
