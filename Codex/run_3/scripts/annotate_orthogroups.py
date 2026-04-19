from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


BLAST_PREFIX_RE = re.compile(r"^(?:sp|tr)\|[^|]+\|[^ ]+\s+")
RELEVANCE_ORDER = {"high": 0, "possible": 1, "low": 2}


def _top_items(counter: Counter[str], limit: int = 5) -> list[dict[str, int | str]]:
    return [
        {"description": description, "count": count}
        for description, count in counter.most_common(limit)
    ]


def _iter_fasta_lengths(path: Path) -> Iterable[int]:
    current: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    yield len("".join(current))
                    current = []
                continue
            current.append(line)
    if current:
        yield len("".join(current))


def parse_blast_file(path: Path) -> dict[str, object]:
    description_counts: Counter[str] = Counter()
    hit_count = 0

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            description = fields[2].split(" OS=")[0]
            description = BLAST_PREFIX_RE.sub("", description)
            description_counts[description] += 1
            hit_count += 1

    return {
        "hit_count": hit_count,
        "top_descriptions": _top_items(description_counts),
    }


def parse_pfam_tbl(path: Path) -> dict[str, object]:
    domain_counts: Counter[str] = Counter()
    domain_hit_count = 0

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 19:
                continue
            description = " ".join(fields[18:])
            domain_counts[description] += 1
            domain_hit_count += 1

    return {
        "domain_hit_count": domain_hit_count,
        "top_domains": _top_items(domain_counts),
    }


def load_signalp_results(directory: Path) -> dict[str, dict[str, object]]:
    grouped: dict[str, dict[str, object]] = defaultdict(
        lambda: {"positive_count": 0, "total_count": 0, "prediction_counts": Counter()}
    )

    for path in sorted(directory.glob("signalp_batch_*.output.json")):
        with path.open() as handle:
            sequences = json.load(handle)["SEQUENCES"]
        for seq_id, record in sequences.items():
            orthogroup = seq_id.split("__", 1)[0]
            prediction = record.get("Prediction", "Other")
            grouped[orthogroup]["total_count"] += 1
            grouped[orthogroup]["prediction_counts"][prediction] += 1
            if prediction != "Other":
                grouped[orthogroup]["positive_count"] += 1

    return dict(grouped)


def load_deeptmhmm_results(directory: Path) -> dict[str, dict[str, object]]:
    grouped: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "total_count": 0,
            "tm_sequence_count": 0,
            "signal_peptide_count": 0,
            "type_counts": Counter(),
        }
    )

    for path in sorted(directory.glob("*.predicted_topologies.3line")):
        lines = path.read_text().splitlines()
        for index in range(0, len(lines), 3):
            if index + 2 >= len(lines):
                break
            header = lines[index]
            topology = lines[index + 2].strip()
            seq_id = header[1:].split()[0]
            orthogroup = seq_id.split("__", 1)[0]
            topology_type = header.split("|")[-1].strip()

            grouped[orthogroup]["total_count"] += 1
            grouped[orthogroup]["type_counts"][topology_type] += 1
            if "M" in topology:
                grouped[orthogroup]["tm_sequence_count"] += 1
            if "S" in topology:
                grouped[orthogroup]["signal_peptide_count"] += 1

    return dict(grouped)


def build_evidence_rows(root: Path) -> list[dict[str, object]]:
    signalp = load_signalp_results(root / "signalp_results")
    deeptmhmm = load_deeptmhmm_results(root / "deeptmhmm_results")
    rows: list[dict[str, object]] = []

    for fasta_path in sorted(root.glob("OG*.fa")):
        orthogroup = fasta_path.stem
        lengths = list(_iter_fasta_lengths(fasta_path))
        blast = parse_blast_file(root / "blast_out" / f"{orthogroup}.blast.tsv")
        pfam = parse_pfam_tbl(root / "hmm_out" / f"{orthogroup}.pfam.tbl")
        signalp_summary = signalp.get(
            orthogroup,
            {"positive_count": 0, "total_count": 0, "prediction_counts": Counter()},
        )
        deeptmhmm_summary = deeptmhmm.get(
            orthogroup,
            {
                "total_count": 0,
                "tm_sequence_count": 0,
                "signal_peptide_count": 0,
                "type_counts": Counter(),
            },
        )

        rows.append(
            {
                "orthogroup": orthogroup,
                "sequence_count": len(lengths),
                "min_length": min(lengths),
                "max_length": max(lengths),
                "blast_hit_count": blast["hit_count"],
                "blast_top_descriptions": blast["top_descriptions"],
                "top_blast_description": (
                    blast["top_descriptions"][0]["description"]
                    if blast["top_descriptions"]
                    else ""
                ),
                "top_blast_count": (
                    blast["top_descriptions"][0]["count"] if blast["top_descriptions"] else 0
                ),
                "pfam_hit_count": pfam["domain_hit_count"],
                "pfam_top_domains": pfam["top_domains"],
                "top_pfam_description": (
                    pfam["top_domains"][0]["description"] if pfam["top_domains"] else ""
                ),
                "top_pfam_count": (
                    pfam["top_domains"][0]["count"] if pfam["top_domains"] else 0
                ),
                "signalp_positive_count": signalp_summary["positive_count"],
                "signalp_total_count": signalp_summary["total_count"],
                "deeptmhmm_total_count": deeptmhmm_summary["total_count"],
                "deeptmhmm_tm_count": deeptmhmm_summary["tm_sequence_count"],
                "deeptmhmm_sp_count": deeptmhmm_summary["signal_peptide_count"],
                "deeptmhmm_type_counts": dict(deeptmhmm_summary["type_counts"]),
            }
        )

    return rows


def load_manual_annotations(path: Path) -> dict[str, dict[str, str]]:
    annotations: dict[str, dict[str, str]] = {}

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            orthogroup = row["orthogroup"]
            annotations[orthogroup] = row

    return annotations


def build_annotated_rows(root: Path, annotation_path: Path) -> list[dict[str, object]]:
    evidence_rows = build_evidence_rows(root)
    annotations = load_manual_annotations(annotation_path)
    evidence_by_og = {row["orthogroup"]: row for row in evidence_rows}

    missing_annotations = sorted(
        orthogroup for orthogroup in evidence_by_og if orthogroup not in annotations
    )
    if missing_annotations:
        raise ValueError(
            "Missing manual annotations for: " + ", ".join(missing_annotations)
        )

    extra_annotations = sorted(
        orthogroup for orthogroup in annotations if orthogroup not in evidence_by_og
    )
    if extra_annotations:
        raise ValueError(
            "Manual annotations do not match input FASTA files: "
            + ", ".join(extra_annotations)
        )

    annotated_rows: list[dict[str, object]] = []
    for orthogroup in sorted(evidence_by_og):
        row = dict(evidence_by_og[orthogroup])
        row.update(annotations[orthogroup])
        row["length_range"] = f"{row['min_length']}-{row['max_length']}"
        row["signalp_fraction"] = (
            f"{row['signalp_positive_count']}/{row['signalp_total_count']}"
        )
        row["deeptmhmm_type_summary"] = _format_type_counts(row["deeptmhmm_type_counts"])
        row["blast_summary"] = _format_item_summary(row["blast_top_descriptions"])
        row["pfam_summary"] = _format_item_summary(row["pfam_top_domains"])
        annotated_rows.append(row)

    annotated_rows.sort(
        key=lambda row: (RELEVANCE_ORDER[row["calcification_relevance"]], row["orthogroup"])
    )
    return annotated_rows


def _format_item_summary(items: list[dict[str, object]]) -> str:
    if not items:
        return "-"
    return "; ".join(f"{item['description']} ({item['count']})" for item in items[:3])


def _format_type_counts(type_counts: dict[str, int]) -> str:
    if not type_counts:
        return "-"
    ordered = sorted(type_counts.items(), key=lambda item: item[0])
    return ", ".join(f"{name}:{count}" for name, count in ordered)


def render_markdown_report(rows: list[dict[str, object]]) -> str:
    high_rows = [row for row in rows if row["calcification_relevance"] == "high"]
    possible_rows = [row for row in rows if row["calcification_relevance"] == "possible"]
    low_rows = [row for row in rows if row["calcification_relevance"] == "low"]

    lines = [
        "# Calcification-Relevant Orthogroup Annotations",
        "",
        (
            "This report summarizes 73 orthogroups using BLAST, PFAM, DeepTMHMM, "
            "and SignalP evidence, then assigns an evidence-based functional call "
            "and an inferred calcification relevance tier."
        ),
        "",
        "## Summary",
        "",
        f"- High relevance: {len(high_rows)} orthogroups",
        f"- Possible/indirect relevance: {len(possible_rows)} orthogroups",
        f"- Low/unclear relevance: {len(low_rows)} orthogroups",
        "",
        "## High-priority calcification candidates",
        "",
    ]

    for row in high_rows:
        lines.append(
            f"- {row['orthogroup']}: {row['predicted_function']}. "
            f"Evidence: BLAST {row['blast_summary']}; PFAM {row['pfam_summary']}; "
            f"DeepTMHMM {row['deeptmhmm_type_summary']}; SignalP {row['signalp_fraction']}. "
            f"{row['calcification_rationale']} "
            f"(confidence: {row['annotation_confidence']})."
        )

    lines.extend(["", "## Possible or indirect candidates", ""])
    for row in possible_rows:
        lines.append(
            f"- {row['orthogroup']}: {row['predicted_function']}. "
            f"Evidence: BLAST {row['blast_summary']}; PFAM {row['pfam_summary']}; "
            f"DeepTMHMM {row['deeptmhmm_type_summary']}; SignalP {row['signalp_fraction']}. "
            f"{row['calcification_rationale']} "
            f"(confidence: {row['annotation_confidence']})."
        )

    lines.extend(
        [
            "",
            "## Low-priority or unresolved candidates",
            "",
            "These orthogroups look more like housekeeping proteins or remain too poorly "
            "resolved for a specific calcification hypothesis:",
            "",
        ]
    )
    for row in low_rows:
        lines.append(
            f"- {row['orthogroup']}: {row['predicted_function']} "
            f"(confidence: {row['annotation_confidence']})."
        )

    lines.extend(
        [
            "",
            "## Full Annotation Table",
            "",
            "| Orthogroup | Predicted function | BLAST summary | PFAM summary | DeepTMHMM | SignalP | Relevance | Confidence |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['orthogroup']} | {row['predicted_function']} | "
            f"{row['blast_summary']} | {row['pfam_summary']} | "
            f"{row['deeptmhmm_type_summary']} | {row['signalp_fraction']} | "
            f"{row['calcification_relevance']} | {row['annotation_confidence']} |"
        )

    lines.append("")
    return "\n".join(lines)


def write_outputs(rows: list[dict[str, object]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "orthogroup_annotation_summary.tsv"
    report_path = output_dir / "orthogroup_annotation_report.md"

    fieldnames = [
        "orthogroup",
        "sequence_count",
        "length_range",
        "top_blast_description",
        "blast_summary",
        "top_pfam_description",
        "pfam_summary",
        "deeptmhmm_type_summary",
        "signalp_fraction",
        "predicted_function",
        "function_basis",
        "calcification_relevance",
        "calcification_rationale",
        "annotation_confidence",
    ]

    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fieldnames})

    report_path.write_text(render_markdown_report(rows))


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    default_root = script_dir.parent
    parser = argparse.ArgumentParser(
        description="Annotate orthogroup functions and calcification relevance."
    )
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument(
        "--annotations",
        type=Path,
        default=default_root / "annotations" / "orthogroup_manual_annotations.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_root / "annotation_results",
    )
    args = parser.parse_args()

    rows = build_annotated_rows(args.root, args.annotations)
    write_outputs(rows, args.output_dir)
    print(f"Wrote {len(rows)} orthogroup annotations to {args.output_dir}")


if __name__ == "__main__":
    main()
