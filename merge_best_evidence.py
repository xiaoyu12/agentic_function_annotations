#!/usr/bin/env python3
"""Merge orthogroup annotation evidence into final tables and a report.

This script is a Python, standard-library implementation of the local
orthogroup annotation merge workflow. It combines:

- local evidence from FASTA, BLAST, Pfam/HMMER, SignalP, and DeepTMHMM;
- nine independent AI annotation runs from the Claude, Claude_code, and Codex
  output folders;
- curated final calls from comparison_curated_final_calls.psv; and
- literature-based updates from Dedman et al. (2024) for coccolith matrix
  proteomics.

The curated final calls remain the authority for biological interpretation.
The deterministic evidence tables make those calls auditable and reproducible.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from statistics import mean, stdev
from typing import Any, Iterable


RUNS = [
    "Claude_app_run1",
    "Claude_app_run2",
    "Claude_app_run3",
    "Claude_code_run1",
    "Claude_code_run2",
    "Claude_code_run3",
    "Codex_app_run1",
    "Codex_app_run2",
    "Codex_app_run3",
]


EVIDENCE_FIELDS = [
    "orthogroup",
    "sequence_count",
    "species_count",
    "length_range",
    "acidic_percent",
    "cys_percent",
    "blast_hit_queries",
    "top_blast",
    "pfam_hit_queries",
    "top_pfam",
    "signalp_positive",
    "signalp_fraction",
    "deeptmhmm_tm_or_beta_count",
    "deeptmhmm_tm_fraction",
    "deeptmhmm_summary",
    "mean_tmrs",
    "max_tmrs",
]


FINAL_FIELDS = [
    "orthogroup",
    "sequence_count",
    "species_count",
    "final_function",
    "final_calcification_relevance",
    "final_rationale",
    "evidence_strength_note",
    "signalp_positive",
    "deeptmhmm_tm_or_beta_count",
    "blast_hit_queries",
    "pfam_hit_queries",
    "top_blast",
    "top_pfam",
    "deeptmhmm_summary",
    "model_consensus_mean_score",
    "model_disagreement_sd",
]


COMPARISON_FIELDS = [
    "orthogroup",
    "final_relevance",
    "final_function",
    "model_score_mean",
    "model_score_sd",
    "model_score_range",
    "overcalled_by_runs",
    "undercalled_by_runs",
    "audit_note",
    *RUNS,
]


METHOD_FIELDS = ["run", "labels", "mean_score", "score_sd", "missing_ogs"]


DEDMP = "Dedman et al. (2024)"


def clean_text(text: Any) -> str:
    if text is None:
        return ""
    value = str(text)
    if not value.strip():
        return ""
    value = value.replace("**", "")
    value = re.sub(r"<[^>]+>", "", value)
    value = re.sub(r"\s+", " ", value)
    return value.strip()


def to_float(value: Any) -> float:
    if value is None:
        return math.nan
    text = str(value).strip()
    if not text or text.upper() == "NA":
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def format_float(value: float, digits: int = 2) -> str:
    if value is None or math.isnan(value):
        return "nan"
    rounded = round(float(value), digits)
    if rounded.is_integer():
        return str(int(rounded))
    return f"{rounded:.{digits}f}".rstrip("0").rstrip(".")


def format_evalue(value: float) -> str:
    if value is None or math.isnan(value):
        return "nan"
    if value == 0:
        return "0"
    if abs(value) < 1e-3 or abs(value) >= 1e4:
        return f"{value:.3g}"
    return f"{value:.6f}".rstrip("0").rstrip(".")


def clean_blast(description: str) -> str:
    text = clean_text(description)
    text = re.sub(r"^sp\|[^|]+\|[^\s]+\s+", "", text)
    text = re.sub(r"^tr\|[^|]+\|[^\s]+\s+", "", text)
    text = re.sub(r"\sOS=.*$", "", text)
    text = re.sub(r"\sOX=\d+.*$", "", text)
    return clean_text(text)


def normalize_relevance(raw: Any) -> tuple[str, float]:
    text = clean_text(raw).lower()
    if not text:
        return "missing", math.nan
    text = re.sub(r"relevance:\s*", "", text)
    text = re.sub(r"[()]", " ", text)
    text = re.sub(r"\s+", " ", text)
    if re.search(r"strong|high priority|high", text):
        return "high", 4.0
    if re.search(r"medium.high|moderate.high", text):
        return "medium-high", 3.5
    if re.search(r"medium|moderate|plausible indirect|possible", text):
        return "moderate/possible", 2.5
    if re.search(r"low.medium|weak|candidate|unknown", text):
        return "watchlist/uncertain", 1.5
    if re.search(r"background|unlikely|low|unclear", text):
        return "low", 0.5
    return text, 1.0


def final_score(label: str) -> float:
    text = clean_text(label).lower()
    if text.startswith("high"):
        return 4.0
    if text.startswith("moderate"):
        return 2.5
    if text.startswith("watchlist"):
        return 1.5
    if text.startswith("low"):
        return 0.5
    return 1.0


def score_stats(scores: Iterable[float]) -> dict[str, float]:
    valid = [score for score in scores if not math.isnan(score)]
    if not valid:
        return {"mean": math.nan, "sd": math.nan, "min": math.nan, "max": math.nan}
    return {
        "mean": round(mean(valid), 2),
        "sd": round(stdev(valid), 2) if len(valid) > 1 else 0.0,
        "min": min(valid),
        "max": max(valid),
    }


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8", errors="replace").splitlines()


def get_og(row: dict[str, str]) -> str:
    for field in ("OG", "og", "orthogroup", "Orthogroup"):
        if field in row:
            return str(row[field])
    return ""


def first_field(row: dict[str, str], fields: list[str]) -> str:
    for field in fields:
        if field in row and row[field] is not None:
            return str(row[field])
    return ""


def add_claim(
    claims: dict[str, dict[str, dict[str, Any]]],
    run: str,
    og: str,
    relevance: str,
    function: str,
    evidence: str,
) -> None:
    if not og:
        return
    label, score = normalize_relevance(relevance)
    claims.setdefault(og, {})[run] = {
        "raw": clean_text(relevance),
        "label": label,
        "score": score,
        "function": clean_text(function),
        "evidence": clean_text(evidence),
    }


def import_tsv_claims(
    claims: dict[str, dict[str, dict[str, Any]]],
    run: str,
    path: Path,
    relevance_fields: list[str],
    function_fields: list[str],
    evidence_fields: list[str],
) -> None:
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            add_claim(
                claims,
                run,
                get_og(row),
                first_field(row, relevance_fields),
                first_field(row, function_fields),
                first_field(row, evidence_fields),
            )


def import_claude_run1(claims: dict[str, dict[str, dict[str, Any]]], path: Path) -> None:
    og: str | None = None
    section: list[str] = []

    def flush() -> None:
        if not og:
            return
        text = " ".join(section)
        relevance = ""
        match = re.search(r"Relevance:\s*\*\*([^*]+)\*\*", text)
        if match:
            relevance = match.group(1)
        else:
            match = re.search(r"Relevance:\s*([^\.]+)", text)
            if match:
                relevance = match.group(1)
        function = ""
        match = re.search(r"Interpretation:\*\*\s*([^\.]+)", text)
        if match:
            function = match.group(1)
        add_claim(claims, "Claude_app_run1", og, relevance, function, text)

    for line in read_lines(path):
        match = re.match(r"^###\s+(OG\d+)", line)
        if match:
            flush()
            og = match.group(1)
            section = []
        elif og:
            section.append(line)
    flush()


def import_claude_run2(claims: dict[str, dict[str, dict[str, Any]]], path: Path) -> None:
    relevance = ""
    for line in read_lines(path):
        if re.match(r"^##\s+High calcification relevance", line):
            relevance = "high"
            continue
        if re.match(r"^##\s+Medium calcification relevance", line):
            relevance = "medium"
            continue
        if re.match(r"^##\s+Low calcification relevance", line):
            relevance = "low"
            continue
        if re.match(r"^##\s+Background calcification relevance", line):
            relevance = "background"
            continue
        match = re.match(r"^\|\s*(OG\d+)\s*\|", line)
        if match:
            parts = line.strip("|").split("|")
            if len(parts) > 7:
                function = f"{clean_text(parts[6])} / {clean_text(parts[7])}"
                add_claim(claims, "Claude_app_run2", match.group(1), relevance, function, line)


def import_claude_run3(claims: dict[str, dict[str, dict[str, Any]]], path: Path) -> None:
    for line in read_lines(path):
        if re.match(r"^\|\s*-", line):
            continue
        match = re.match(r"^\|\s*\**(OG\d+)\**\s*\|", line)
        if not match:
            continue
        parts = line.strip("|").split("|")
        if len(parts) >= 8:
            og = parts[0].replace("*", "").strip()
            function = f"{clean_text(parts[5])} / {clean_text(parts[6])}"
            add_claim(claims, "Claude_app_run3", og, clean_text(parts[7]), function, line)


def load_claims(root: Path) -> tuple[dict[str, dict[str, dict[str, Any]]], list[str]]:
    claims: dict[str, dict[str, dict[str, Any]]] = {}
    warnings: list[str] = []

    custom_imports = [
        (import_claude_run1, root / "Claude" / "run_1" / "OG_annotation_report.md"),
        (import_claude_run2, root / "Claude" / "run_2" / "annotations" / "og_summary.md"),
        (import_claude_run3, root / "Claude" / "run_3" / "annotation" / "OG_annotation_report.md"),
    ]
    for importer, path in custom_imports:
        if path.exists():
            importer(claims, path)
        else:
            warnings.append(f"Missing model output: {path}")

    tsv_imports = [
        (
            "Claude_code_run1",
            root / "Claude_code" / "run_1" / "annotation_output" / "og_annotations.tsv",
            ["calcification_relevance"],
            ["function_short", "function_long"],
            ["calcification_rationale", "top_pfam", "top_blast"],
        ),
        (
            "Claude_code_run2",
            root / "Claude_code" / "run_2" / "annotation_calcification.tsv",
            ["calcification_relevance"],
            ["primary_annotation"],
            ["rationale", "key_domains_or_hits"],
        ),
        (
            "Claude_code_run3",
            root / "Claude_code" / "run_3" / "orthogroup_annotation.tsv",
            ["calcification_relevance"],
            ["function_summary"],
            ["calcification_rationale", "top_pfam_domains", "top_blast_hits"],
        ),
        (
            "Codex_app_run1",
            root / "Codex" / "run_1" / "annotation_results" / "orthogroup_annotation_summary.tsv",
            ["calcification_relevance"],
            ["inferred_function"],
            ["calcification_rationale", "top_pfam_domains", "top_blast_consensus"],
        ),
        (
            "Codex_app_run2",
            root / "Codex" / "run_2" / "orthogroup_annotation_evidence.tsv",
            ["calcification_relevance"],
            ["function_call"],
            ["calcification_rationale", "top_pfam_summary", "top_blast_summary"],
        ),
        (
            "Codex_app_run3",
            root / "Codex" / "run_3" / "annotations" / "orthogroup_manual_annotations.tsv",
            ["calcification_relevance"],
            ["predicted_function"],
            ["calcification_rationale", "function_basis"],
        ),
    ]
    for run, path, relevance_fields, function_fields, evidence_fields in tsv_imports:
        if path.exists():
            import_tsv_claims(claims, run, path, relevance_fields, function_fields, evidence_fields)
        else:
            warnings.append(f"Missing model output: {path}")

    return claims, warnings


def get_fasta_stats(input_dir: Path) -> dict[str, dict[str, Any]]:
    stats: dict[str, dict[str, Any]] = {}
    for file_path in sorted(input_dir.glob("OG*.fa")):
        og = file_path.stem
        sequences: list[str] = []
        species: Counter[str] = Counter()
        header = ""
        chunks: list[str] = []
        for line in read_lines(file_path):
            if line.startswith(">"):
                if header:
                    sequences.append(re.sub(r"X|\*", "", re.sub(r"[^A-Za-z]", "", "".join(chunks))))
                    parts = header.lstrip(">").split("|")
                    if len(parts) >= 2:
                        species[parts[1]] += 1
                header = line
                chunks = []
            else:
                chunks.append(line.strip())
        if header:
            sequences.append(re.sub(r"X|\*", "", re.sub(r"[^A-Za-z]", "", "".join(chunks))))
            parts = header.lstrip(">").split("|")
            if len(parts) >= 2:
                species[parts[1]] += 1

        lengths = sorted(len(seq) for seq in sequences)
        all_sequence = "".join(sequences)
        stats[og] = {
            "sequence_count": len(sequences),
            "species_count": len(species),
            "mean_length": round(mean(lengths), 1) if lengths else 0,
            "min_length": lengths[0] if lengths else 0,
            "max_length": lengths[-1] if lengths else 0,
            "acidic_percent": round((len(re.findall(r"[DE]", all_sequence)) / len(all_sequence)) * 100, 1)
            if all_sequence
            else 0,
            "cys_percent": round((all_sequence.count("C") / len(all_sequence)) * 100, 1) if all_sequence else 0,
        }
    return stats


def get_blast_evidence(blast_dir: Path) -> dict[str, dict[str, Any]]:
    evidence: dict[str, dict[str, Any]] = {}
    for file_path in blast_dir.glob("OG*.blast.tsv"):
        og = file_path.name.replace(".blast.tsv", "")
        best: dict[str, dict[str, Any]] = {}
        for line in read_lines(file_path):
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            evalue = to_float(parts[5])
            if math.isnan(evalue) or evalue > 1e-3:
                continue
            bits = to_float(parts[6])
            record = {
                "query": parts[0],
                "desc": clean_blast(parts[2]),
                "evalue": evalue,
                "bits": bits,
            }
            old = best.get(parts[0])
            if old is None or evalue < old["evalue"] or (evalue == old["evalue"] and bits > old["bits"]):
                best[parts[0]] = record
        groups = Counter(record["desc"] for record in best.values())
        top = "; ".join(f"{name} ({count})" for name, count in sorted(groups.items(), key=lambda x: (-x[1], x[0]))[:5])
        evidence[og] = {"blast_hit_queries": len(best), "top_blast": top}
    return evidence


def get_pfam_evidence(pfam_dir: Path) -> dict[str, dict[str, Any]]:
    evidence: dict[str, dict[str, Any]] = {}
    for file_path in pfam_dir.glob("OG*.pfam.tbl"):
        og = file_path.name.replace(".pfam.tbl", "")
        records: list[dict[str, Any]] = []
        for line in read_lines(file_path):
            if not line.strip() or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line.strip(), maxsplit=18)
            if len(parts) < 18:
                continue
            evalue = to_float(parts[4])
            if math.isnan(evalue) or evalue > 1e-3:
                continue
            records.append(
                {
                    "domain": parts[0],
                    "query": parts[2],
                    "evalue": evalue,
                    "desc": clean_text(parts[18]) if len(parts) >= 19 else "",
                }
            )
        grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for record in records:
            grouped[record["domain"]].append(record)
        top_parts: list[str] = []
        for domain, group in sorted(grouped.items(), key=lambda x: (-len(x[1]), x[0]))[:6]:
            best = sorted(group, key=lambda x: x["evalue"])[0]
            top_parts.append(f"{domain} ({len(group)}; {best['desc']}; bestE={format_evalue(best['evalue'])})")
        evidence[og] = {
            "pfam_hit_queries": len({record["query"] for record in records}),
            "top_pfam": "; ".join(top_parts),
        }
    return evidence


def get_signalp_evidence(signalp_dir: Path) -> dict[str, dict[str, Any]]:
    counts: dict[str, dict[str, int]] = defaultdict(lambda: {"pos": 0, "total": 0})
    for file_path in signalp_dir.rglob("prediction_results.txt"):
        for line in read_lines(file_path):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            match = re.match(r"^(OG\d+)__", parts[0])
            if not match:
                continue
            og = match.group(1)
            counts[og]["total"] += 1
            if parts[1] != "OTHER":
                counts[og]["pos"] += 1
    return {
        og: {
            "signalp_positive": values["pos"],
            "signalp_total": values["total"],
            "signalp_fraction": round(values["pos"] / values["total"], 3) if values["total"] else 0,
        }
        for og, values in counts.items()
    }


def get_tmhmm_evidence(tmhmm_dir: Path) -> dict[str, dict[str, Any]]:
    labels: dict[str, Counter[str]] = defaultdict(Counter)
    for file_path in tmhmm_dir.glob("*.predicted_topologies.3line"):
        for line in read_lines(file_path):
            match = re.match(r"^>(OG\d+)__\d+\s+\|\s+(\S+)", line)
            if match:
                labels[match.group(1)][match.group(2)] += 1

    tmrs: dict[str, list[int]] = defaultdict(list)
    for file_path in tmhmm_dir.glob("*.TMRs.gff3"):
        for line in read_lines(file_path):
            match = re.match(r"^#\s+(OG\d+__\d+)\s+Number of predicted TMRs:\s+(\d+)", line)
            if match:
                og = re.sub(r"__\d+$", "", match.group(1))
                tmrs[og].append(int(match.group(2)))

    evidence: dict[str, dict[str, Any]] = {}
    for og, label_counts in labels.items():
        total = sum(label_counts.values())
        tm_count = sum(label_counts.get(label, 0) for label in ("TM", "SP+TM", "BETA"))
        summary = "; ".join(f"{label}:{label_counts[label]}" for label in sorted(label_counts))
        og_tmrs = tmrs.get(og, [])
        evidence[og] = {
            "deeptmhmm_summary": summary,
            "tm_or_beta_count": tm_count,
            "tm_fraction": round(tm_count / total, 3) if total else 0,
            "mean_tmrs": round(mean(og_tmrs), 2) if og_tmrs else 0,
            "max_tmrs": max(og_tmrs) if og_tmrs else 0,
        }
    return evidence


def recompute_local_evidence(input_dir: Path) -> dict[str, dict[str, Any]]:
    fasta = get_fasta_stats(input_dir)
    blast = get_blast_evidence(input_dir / "blast_out")
    pfam = get_pfam_evidence(input_dir / "hmm_out")
    signalp = get_signalp_evidence(input_dir / "signalp_results")
    tmhmm = get_tmhmm_evidence(input_dir / "deeptmhmm_results")

    rows: dict[str, dict[str, Any]] = {}
    for og in sorted(fasta):
        fs = fasta[og]
        be = blast.get(og, {})
        pe = pfam.get(og, {})
        se = signalp.get(og, {})
        te = tmhmm.get(og, {})
        rows[og] = {
            "orthogroup": og,
            "sequence_count": fs["sequence_count"],
            "species_count": fs["species_count"],
            "length_range": f"{fs['min_length']}-{fs['max_length']}",
            "acidic_percent": fs["acidic_percent"],
            "cys_percent": fs["cys_percent"],
            "blast_hit_queries": be.get("blast_hit_queries", 0),
            "top_blast": be.get("top_blast", ""),
            "pfam_hit_queries": pe.get("pfam_hit_queries", 0),
            "top_pfam": pe.get("top_pfam", ""),
            "signalp_positive": se.get("signalp_positive", 0),
            "signalp_fraction": se.get("signalp_fraction", 0),
            "deeptmhmm_tm_or_beta_count": te.get("tm_or_beta_count", 0),
            "deeptmhmm_tm_fraction": te.get("tm_fraction", 0),
            "deeptmhmm_summary": te.get("deeptmhmm_summary", ""),
            "mean_tmrs": te.get("mean_tmrs", 0),
            "max_tmrs": te.get("max_tmrs", 0),
        }
    return rows


def apply_literature_updates(final_calls: dict[str, dict[str, str]]) -> None:
    updates = {
        "OG0001976": {
            "final_relevance": "Watchlist",
            "final_rationale": (
                "Strong pentapeptide-repeat support with limited secretion; Dedman et al. (2024) provides "
                "external proteomic support for pentapeptide repeats as recurring coccolith-matrix features, "
                "but this remains a motif-level follow-up rather than a proven Ca/carbonate-binding function."
            ),
            "audit_note": (
                "Literature support justifies watchlist status; Ca/carbonate-binding claims remain overreach "
                "without orthogroup-specific validation."
            ),
        },
        "OG0010867": {
            "final_relevance": "Watchlist",
            "final_rationale": (
                "Pentapeptide-repeat support with weak secretion; Dedman et al. (2024) identifies "
                "pentapeptide-repeat proteins in coccolith matrices across species, supporting follow-up "
                "relevance but not a direct mineral-binding annotation."
            ),
            "audit_note": "Literature support is motif-level and does not prove this orthogroup is a coccolith protein.",
        },
        "OG0022524": {
            "final_relevance": "Watchlist",
            "final_rationale": (
                "Strong pentapeptide-repeat support; Dedman et al. (2024) supports pentapeptide repeats as "
                "coccolith-matrix-associated motifs, but this orthogroup lacks secretion/topology support and "
                "remains a cautious follow-up candidate."
            ),
            "audit_note": "Upgraded only to watchlist because no SignalP or DeepTMHMM secretion support is present.",
        },
        "OG0021347": {
            "final_rationale": (
                "Strong trypsin PFAM plus secretion; Dedman et al. (2024) reports peptidases/protease-regulation "
                "as coccolith-matrix themes, so matrix remodeling is plausible but still indirect."
            ),
        },
        "OG0023594": {
            "final_rationale": (
                "GT8 support is consistent with extracellular polysaccharide biosynthesis, and Dedman et al. "
                "(2024) reports glycosyltransferase/carbohydrate-metabolism proteins in coccolith preparations; "
                "still indirect without orthogroup-specific localization."
            ),
        },
        "OG0024846": {
            "final_rationale": (
                "MUCI70/TOD1 glycosyltransferase-like domain suggests matrix/cell-wall polysaccharide modification; "
                "Dedman et al. (2024) supports carbohydrate-modifying proteins as plausible coccolithogenesis candidates."
            ),
        },
    }
    for og, fields in updates.items():
        if og in final_calls:
            final_calls[og].update(fields)


def load_final_calls(path: Path, apply_literature: bool = True) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        rows = {
            row["orthogroup"]: {
                "orthogroup": row.get("orthogroup", ""),
                "final_function": row.get("final_function", ""),
                "final_relevance": row.get("final_relevance", ""),
                "final_rationale": row.get("final_rationale", ""),
                "audit_note": row.get("audit_note", ""),
            }
            for row in csv.DictReader(handle, delimiter="|")
            if row.get("orthogroup")
        }
    if apply_literature:
        apply_literature_updates(rows)
    return rows


def write_delimited(path: Path, rows: list[dict[str, Any]], fields: list[str], delimiter: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_report(
    final_rows: list[dict[str, Any]],
    method_rows: list[dict[str, Any]],
    evidence_rows: list[dict[str, Any]],
    input_counts: dict[str, Any],
) -> str:
    distribution = Counter(row["final_calcification_relevance"] for row in final_rows)
    distribution_text = "; ".join(f"{label}: {distribution[label]}" for label in sorted(distribution))
    high = ", ".join(row["orthogroup"] for row in final_rows if row["final_calcification_relevance"] == "High")
    moderate = ", ".join(row["orthogroup"] for row in final_rows if row["final_calcification_relevance"] == "Moderate")
    watchlist = ", ".join(row["orthogroup"] for row in final_rows if row["final_calcification_relevance"] == "Watchlist")

    lines = [
        "# Comparison of AI Orthogroup Annotation Outputs",
        "",
        f"Generated on {datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %z')} from the local shared evidence and nine model runs.",
        "",
        "## Inputs Audited",
        "",
        f"- FASTA orthogroups: {input_counts['fasta_orthogroups']}",
        f"- Total protein sequences: {input_counts['total_sequences']}",
        f"- BLAST tables: {input_counts['blast_tables']}",
        f"- PFAM/HMMER tables: {input_counts['pfam_tables']}",
        f"- DeepTMHMM topology batches: {input_counts['tmhmm_batches']}",
        f"- SignalP result batches: {input_counts['signalp_batches']}",
        "",
        "## Final Merged Relevance Distribution",
        "",
        distribution_text,
        "",
        f"High-confidence final candidates: {high}",
        "",
        f"Moderate candidates: {moderate}",
        "",
        f"Watchlist candidates: {watchlist}",
        "",
        "## Method-Level Consistency",
        "",
        "| Run | Normalized label distribution | Mean score | Missing OGs |",
        "| --- | --- | ---: | ---: |",
    ]
    for row in method_rows:
        lines.append(f"| {row['run']} | {row['labels']} | {row['mean_score']} | {row['missing_ogs']} |")
    lines.extend(
        [
            "",
            "The largest source of inconsistency was biological interpretation, not parsing. Keyword/topology-heavy runs promoted many secreted or membrane proteins without calcification-specific domains, while conservative runs left those as low or watchlist candidates.",
            "",
            "## Literature Review Update: Dedman et al. Coccolith Matrix Proteomics",
            "",
            'The attached paper, Dedman et al. (2024), Scientific Reports, "Exploring proteins within the coccolith matrix" (DOI: 10.1038/s41598-024-83052-9), adds useful external proteomic context. The study analyzed cleaned coccolith material from Gephyrocapsa huxleyi, Gephyrocapsa oceanica, and Coccolithus braarudii, and compared conserved protein features with shell or skeletal matrix proteins from other marine calcifiers.',
            "",
            "The strongest annotation impact is on pentapeptide-repeat orthogroups. Dedman et al. report pentapeptide-repeat proteins in all three examined coccolithophore species and explicitly propose this repetitive structural motif as a coccolith-matrix-associated feature worthy of targeted functional testing. Therefore OG0001976, OG0010867, and OG0022524 have been moved from Low to Watchlist. They were not promoted to Moderate or High because the support is motif-level rather than orthogroup-specific, secretion support is weak or absent in the local data, and the paper itself treats the mechanism as unresolved.",
            "",
            "The paper also supports several existing cautious annotations without changing their tiers. Its discussion of carbohydrate-modifying proteins and CAP/baseplate chemistry reinforces the Moderate tier for glycosyltransferase-like candidates such as OG0023594 and OG0024846. Its report of coccolith peptidases and protease-regulation features supports the plausibility of OG0021347 as a secreted matrix-remodeling candidate, while retaining the Moderate tier because trypsin-family enzymes are not calcification-specific. The paper does not justify upgrading generic housekeeping, SMC, histone, ribosomal, or carbon-metabolism annotations, and it cautions that low-abundance matrix proteins require localization or genetic validation.",
            "",
            "## Biological Correctness And Overclaim Patterns",
            "",
            "1. Sulfation and sulfated-glycan biology was the most reproducible direct signal. OG0017138 and OG0020703 remain high-confidence final candidates.",
            "2. Glycosyltransferase-like OGs are plausible but mostly indirect. OG0009816, OG0011061, OG0017305, OG0023594, and OG0024846 are retained as moderate; OG0011061 and OG0017305 were downgraded from some high calls because domain coverage is limited.",
            "3. Secreted and membrane orphan families are important follow-up targets but should not be annotated as known calcification proteins. They are separated into moderate or watchlist tiers depending on topology strength.",
            "4. Pentapeptide-repeat OGs were often over-interpreted in the original runs. The Dedman et al. proteomics paper adds external support for pentapeptide repeats as coccolith-matrix-associated motifs, so those OGs are now watchlist candidates, but direct Ca/carbonate-binding or proven calcification-function claims remain overreach.",
            "5. BLAST hits to giant/repetitive proteins caused misleading functions in multiple outputs. RPB1 for OG0009301, Piccolo for OG0011061/OG0015153, and collagen calls based on one member of OG0016203 need domain/topology corroboration.",
            "",
            "## Output Files",
            "",
            "- `recomputed_evidence.tsv`: direct re-aggregation of FASTA, BLAST, PFAM, SignalP, and DeepTMHMM.",
            "- `model_run_comparison.tsv`: one row per OG with all nine relevance calls, disagreement metrics, and audit notes.",
            "- `method_consistency_summary.tsv`: normalized tier distribution for each model run.",
            "- `final_orthogroup_annotations.tsv`: final merged annotation table.",
            "",
            "## Final Annotation Table",
            "",
            "| OG | n | Final function | Relevance | Rationale |",
            "| --- | ---: | --- | --- | --- |",
        ]
    )
    for row in final_rows:
        function = clean_text(row["final_function"]).replace("|", "/")
        rationale = clean_text(row["final_rationale"]).replace("|", "/")
        lines.append(
            f"| {row['orthogroup']} | {row['sequence_count']} | {function} | {row['final_calcification_relevance']} | {rationale} |"
        )
    return "\n".join(lines) + "\n"


def build_outputs(
    root: Path,
    input_dir: Path,
    curated_calls_path: Path,
    out_dir: Path,
    apply_literature: bool = True,
) -> dict[str, Any]:
    claims, warnings = load_claims(root)
    final_calls = load_final_calls(curated_calls_path, apply_literature=apply_literature)
    evidence_by_og = recompute_local_evidence(input_dir)
    ogs = sorted(evidence_by_og)

    missing_final_calls = [og for og in ogs if og not in final_calls]
    if missing_final_calls:
        missing = ", ".join(missing_final_calls[:10])
        raise ValueError(f"Missing curated final calls for {len(missing_final_calls)} orthogroups: {missing}")

    evidence_rows: list[dict[str, Any]] = [evidence_by_og[og] for og in ogs]
    comparison_rows: list[dict[str, Any]] = []
    final_rows: list[dict[str, Any]] = []

    for og in ogs:
        ev = evidence_by_og[og]
        final = final_calls[og]
        scores = [claims[og][run]["score"] for run in RUNS if og in claims and run in claims[og]]
        stats = score_stats(scores)
        fscore = final_score(final["final_relevance"])
        overcalled: list[str] = []
        undercalled: list[str] = []
        for run in RUNS:
            if og not in claims or run not in claims[og]:
                continue
            score = claims[og][run]["score"]
            if math.isnan(score):
                continue
            if score >= fscore + 1.0:
                overcalled.append(run)
            if score <= fscore - 1.0:
                undercalled.append(run)

        comparison = {
            "orthogroup": og,
            "final_relevance": final["final_relevance"],
            "final_function": final["final_function"],
            "model_score_mean": format_float(stats["mean"]),
            "model_score_sd": format_float(stats["sd"]),
            "model_score_range": f"{format_float(stats['min'])}-{format_float(stats['max'])}",
            "overcalled_by_runs": ";".join(overcalled),
            "undercalled_by_runs": ";".join(undercalled),
            "audit_note": final["audit_note"],
        }
        for run in RUNS:
            comparison[run] = claims.get(og, {}).get(run, {}).get("raw", "")
        comparison_rows.append(comparison)

        final_rows.append(
            {
                "orthogroup": og,
                "sequence_count": ev["sequence_count"],
                "species_count": ev["species_count"],
                "final_function": final["final_function"],
                "final_calcification_relevance": final["final_relevance"],
                "final_rationale": final["final_rationale"],
                "evidence_strength_note": final["audit_note"],
                "signalp_positive": ev["signalp_positive"],
                "deeptmhmm_tm_or_beta_count": ev["deeptmhmm_tm_or_beta_count"],
                "blast_hit_queries": ev["blast_hit_queries"],
                "pfam_hit_queries": ev["pfam_hit_queries"],
                "top_blast": ev["top_blast"],
                "top_pfam": ev["top_pfam"],
                "deeptmhmm_summary": ev["deeptmhmm_summary"],
                "model_consensus_mean_score": format_float(stats["mean"]),
                "model_disagreement_sd": format_float(stats["sd"]),
            }
        )

    method_rows: list[dict[str, Any]] = []
    for run in RUNS:
        labels: list[str] = []
        scores: list[float] = []
        for og in ogs:
            claim = claims.get(og, {}).get(run)
            if claim is None:
                labels.append("missing")
            else:
                labels.append(claim["label"])
                scores.append(claim["score"])
        stats = score_stats(scores)
        label_counts = Counter(labels)
        method_rows.append(
            {
                "run": run,
                "labels": "; ".join(f"{label}={label_counts[label]}" for label in sorted(label_counts)),
                "mean_score": format_float(stats["mean"]),
                "score_sd": format_float(stats["sd"]),
                "missing_ogs": label_counts.get("missing", 0),
            }
        )

    input_counts = {
        "fasta_orthogroups": len(ogs),
        "total_sequences": sum(row["sequence_count"] for row in evidence_rows),
        "blast_tables": len(list((input_dir / "blast_out").glob("OG*.blast.tsv"))),
        "pfam_tables": len(list((input_dir / "hmm_out").glob("OG*.pfam.tbl"))),
        "tmhmm_batches": len(list((input_dir / "deeptmhmm_results").glob("*.predicted_topologies.3line"))),
        "signalp_batches": len(list((input_dir / "signalp_results").rglob("prediction_results.txt"))),
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    write_delimited(out_dir / "recomputed_evidence.tsv", evidence_rows, EVIDENCE_FIELDS)
    write_delimited(out_dir / "model_run_comparison.tsv", comparison_rows, COMPARISON_FIELDS)
    write_delimited(out_dir / "method_consistency_summary.tsv", method_rows, METHOD_FIELDS)
    write_delimited(out_dir / "final_orthogroup_annotations.tsv", final_rows, FINAL_FIELDS)
    (out_dir / "comparison_and_final_annotation_report.md").write_text(
        build_report(final_rows, method_rows, evidence_rows, input_counts),
        encoding="utf-8",
    )

    return {
        "orthogroups": len(ogs),
        "total_sequences": input_counts["total_sequences"],
        "distribution": dict(Counter(row["final_calcification_relevance"] for row in final_rows)),
        "warnings": warnings,
        "out_dir": str(out_dir),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="Project/output root directory.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Input evidence directory. Defaults to ROOT/Orthogroups.calcifying_loose_fastas.",
    )
    parser.add_argument(
        "--curated-calls",
        type=Path,
        default=None,
        help="Curated final calls PSV. Defaults to ROOT/comparison_curated_final_calls.psv.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to ROOT/comparison_merged_annotation.",
    )
    parser.add_argument(
        "--no-literature-updates",
        action="store_true",
        help="Do not apply the Dedman et al. literature update layer to curated calls.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    input_dir = (args.input_dir or root / "Orthogroups.calcifying_loose_fastas").resolve()
    curated_calls = (args.curated_calls or root / "comparison_curated_final_calls.psv").resolve()
    out_dir = (args.out_dir or root / "comparison_merged_annotation").resolve()

    result = build_outputs(
        root=root,
        input_dir=input_dir,
        curated_calls_path=curated_calls,
        out_dir=out_dir,
        apply_literature=not args.no_literature_updates,
    )

    distribution = "; ".join(f"{key}={value}" for key, value in sorted(result["distribution"].items()))
    print(f"Wrote comparison outputs to {result['out_dir']}")
    print(f"Orthogroups: {result['orthogroups']}; sequences: {result['total_sequences']}; {distribution}")
    for warning in result["warnings"]:
        print(f"WARNING: {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
