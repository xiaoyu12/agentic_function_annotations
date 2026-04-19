from __future__ import annotations

import csv
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


DIRECT_CALCIFICATION_KEYWORDS = (
    "carbonic anhydrase",
    "bicarbonate",
    "carbonate",
    "calcium",
    "cation channel",
    "calmodulin",
    "proton pump",
    "proton transporter",
    "v-type atpase",
    "caax amino terminal protease",
    "anion transporter",
    "transporter",
    "channel",
    "exchanger",
    "antiporter",
    "symporter",
    "atpase",
    "ion homeostasis",
)

PLAUSIBLE_CALCIFICATION_KEYWORDS = (
    "secreted",
    "signal peptide",
    "extracellular",
    "cell surface",
    "golgi",
    "vesicle",
    "trafficking",
    "snare",
    "rab",
    "syntaxin",
    "coatomer",
    "glycosyltransferase",
    "glycosidase",
    "kinase",
    "receptor",
    "adhesion",
    "dynein",
    "kinesin",
    "microtubule",
    "actin",
    "cytoskeleton",
    "coiled-coil",
    "redox",
    "oxidase",
    "peroxidase",
    "glycosyltransferase",
    "mannosyltransferase",
    "sulfatase",
    "sulfotransferase",
    "protease",
    "ubiquitin",
    "phosphatase",
)

LOW_RELEVANCE_KEYWORDS = (
    "ribosomal",
    "mitotic checkpoint",
    "histone",
    "endonuclease",
    "dna repair",
    "rna",
    "translation",
    "splicing",
    "proteasome",
    "tubulin-folding",
    "replication",
)

GENERIC_DESCRIPTION_PATTERNS = (
    "uncharacterized protein",
    "hypothetical protein",
    "predicted protein",
    "unknown protein",
    "domain of unknown function",
)

GENERIC_PFAM_PATTERNS = (
    "aaa domain",
    "fad binding domain",
    "oxidoreductase nad-binding domain",
    "viral superfamily 1",
    "short c-terminal domain",
    "pdz domain",
    "beta-propeller",
)


def orthogroup_from_surrogate_id(sequence_id: str) -> str:
    return sequence_id.split("__", 1)[0]


def normalize_blast_description(description: str) -> str:
    text = description.strip()
    if text.startswith(("sp|", "tr|", "gi|")) and " " in text:
        text = text.split(" ", 1)[1]
    text = re.split(r"\sOS=", text, maxsplit=1)[0]
    text = re.sub(r"\s+\(Fragment\)$", "", text, flags=re.IGNORECASE)
    text = re.sub(r"\s+(putative|probable|predicted)$", "", text, flags=re.IGNORECASE)
    return re.sub(r"\s+", " ", text).strip()


def parse_fasta(path: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    header: str | None = None
    sequence_chunks: list[str] = []
    orthogroup = path.stem
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(sequence_chunks)
                    records.append(
                        {
                            "orthogroup": orthogroup,
                            "original_id": header,
                            "sequence": sequence,
                            "length": len(sequence),
                            "surrogate_id": f"{orthogroup}__{len(records) + 1:06d}",
                        }
                    )
                header = line[1:].strip().split()[0]
                sequence_chunks = []
            else:
                sequence_chunks.append(line.strip())
    if header is not None:
        sequence = "".join(sequence_chunks)
        records.append(
            {
                "orthogroup": orthogroup,
                "original_id": header,
                "sequence": sequence,
                "length": len(sequence),
                "surrogate_id": f"{orthogroup}__{len(records) + 1:06d}",
            }
        )
    return records


def parse_pfam_tbl(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split(maxsplit=18)
            if len(fields) < 19:
                continue
            rows.append(
                {
                    "domain": fields[0],
                    "accession": fields[1],
                    "query": fields[2],
                    "full_evalue": float(fields[4]),
                    "full_score": float(fields[5]),
                    "domain_evalue": float(fields[7]),
                    "description": fields[18],
                }
            )
    return rows


def parse_deeptmhmm_3line(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            topology = handle.readline()
            if not sequence or not topology:
                break
            header = header.strip()
            if not header:
                continue
            header_body = header[1:] if header.startswith(">") else header
            sequence_id, label = [part.strip() for part in header_body.split("|", maxsplit=1)]
            rows.append(
                {
                    "sequence_id": sequence_id,
                    "label": label,
                    "length": len(sequence.strip()),
                    "topology": topology.strip(),
                }
            )
    return rows


def parse_signalp_json(path: Path) -> list[dict[str, object]]:
    payload = json.loads(path.read_text())
    rows: list[dict[str, object]] = []
    for sequence_id, record in payload["SEQUENCES"].items():
        rows.append(
            {
                "sequence_id": sequence_id,
                "prediction": record["Prediction"],
                "cleavage_site": record["CS_pos"],
            }
        )
    return rows


def parse_blast_file(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 7:
                continue
            description = fields[2]
            rows.append(
                {
                    "query": fields[0],
                    "subject_id": fields[1],
                    "description": description,
                    "normalized_description": normalize_blast_description(description),
                    "pident": float(fields[3]),
                    "alignment_length": int(fields[4]),
                    "evalue": float(fields[5]),
                    "bitscore": float(fields[6]),
                }
            )
    return rows


def summarize_counter(counter: Counter[str], limit: int = 5) -> str:
    parts: list[str] = []
    for label, count in counter.most_common(limit):
        parts.append(f"{label} ({count})")
    return "; ".join(parts)


def median_or_zero(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.median(values) if values else 0.0


def is_generic_description(text: str) -> bool:
    lowered = text.lower()
    return any(pattern in lowered for pattern in GENERIC_DESCRIPTION_PATTERNS)


def clean_pfam_description(description: str) -> str:
    text = re.sub(r"\s+", " ", description).strip()
    if text.lower().endswith("-like domain"):
        return text[:-12] + "-like protein"
    if text.lower().endswith(" domain"):
        return f"protein containing {text}"
    return text


def top_informative_label(counter: Counter[str]) -> str:
    for label, _ in counter.most_common():
        if not is_generic_description(label):
            return label
    return counter.most_common(1)[0][0] if counter else ""


def is_generic_pfam_description(description: str) -> bool:
    lowered = description.lower()
    return any(pattern in lowered for pattern in GENERIC_PFAM_PATTERNS)


def is_membrane_compatible_label(label: str) -> bool:
    lowered = label.lower()
    return any(
        token in lowered
        for token in (
            "membrane",
            "transmembrane",
            "receptor",
            "channel",
            "transporter",
            "permease",
            "gpcr",
            "sur7",
            "pali",
            "warthog",
            "hedgehog",
            "rhd3",
            "sey1",
        )
    )


def family_label_from_pfam(description: str, signalp_fraction: float, tm_fraction: float) -> str:
    lowered = description.lower()
    if not description:
        return ""
    if "hint module" in lowered or "hint domain" in lowered:
        return "Hedgehog/Warthog-like Hint-domain membrane protein"
    if "fg-gap" in lowered:
        return "FG-GAP repeat-containing receptor-like membrane protein"
    if "ankyrin repeat" in lowered:
        return "ankyrin-repeat protein"
    if "kelch" in lowered:
        return "kelch-repeat beta-propeller protein"
    if "leucine rich repeat" in lowered:
        return "leucine-rich repeat protein"
    if "exostosin gt47" in lowered:
        return "GT47 exostosin-like membrane glycosyltransferase"
    if "glycosyl transferase family 8" in lowered:
        return "glycosyltransferase family 8 protein"
    if "mannosyltransferase" in lowered:
        return "mannosyltransferase-like enzyme"
    if "sulfotransferase" in lowered:
        return "sulfotransferase-like enzyme"
    if "sulphotransferase" in lowered:
        return "sulfotransferase-like enzyme"
    if "sulfatase" in lowered:
        return "sulfatase-family hydrolase"
    if "sur7/pali" in lowered:
        return "SUR7/PalI-like membrane protein"
    if "short c-terminal domain" in lowered:
        return "proline-rich membrane-anchored protein"
    if "pentapeptide" in lowered:
        return "pentapeptide repeat protein"
    if "2og-fe(ii) oxygenase" in lowered:
        return "2OG-Fe(II) oxygenase"
    if "trypsin" in lowered:
        return "trypsin-like serine protease"
    if "amidase" in lowered:
        return "amidase-family enzyme"
    if "oxidoreductase fad-binding domain" in lowered or "nad(p)-binding rossmann-like domain" in lowered:
        return "flavin-dependent oxidoreductase"
    if "cyclic nucleotide-binding domain" in lowered:
        return "cyclic-nucleotide-binding protein"
    if "serine aminopeptidase, s33" in lowered:
        return "serine aminopeptidase"
    if "pif1-like helicase" in lowered:
        return "PIF1-like helicase"
    if "aaa domain" in lowered:
        return "ATP-dependent helicase-like protein"
    if "f-box-like" in lowered:
        return "F-box protein"
    if "hect-domain" in lowered:
        return "HECT-type E3 ubiquitin ligase"
    if "protein phosphatase 2c" in lowered:
        return "PP2C phosphatase"
    if "pdz domain" in lowered:
        return "PDZ-domain protein"
    if "jmjc domain" in lowered:
        return "JmjC-domain oxygenase"
    if "collagen triple helix repeat" in lowered:
        return "collagen-like repeat protein"
    if "tod1/muci70" in lowered:
        return "MUCI70/TOD1-like glycosyltransferase"
    if "50s ribosome-binding gtpase" in lowered:
        return "ribosome-associated GTPase"
    if "glutathione s-transferase" in lowered:
        return "glutathione S-transferase-like protein"
    return clean_pfam_description(description)


def family_label_from_blast(description: str) -> str:
    lowered = description.lower()
    if not description:
        return ""
    if "warthog protein" in lowered or "protein hedgehog" in lowered:
        return "warthog/hedgehog-like membrane protein"
    if "metabotropic glutamate receptor" in lowered:
        return "GPCR-like receptor"
    if "arylsulfatase" in lowered:
        return "arylsulfatase-like enzyme"
    if "heparan-sulfate 6-o-sulfotransferase" in lowered:
        return "heparan-sulfate sulfotransferase-like enzyme"
    if "protein sey1 homolog" in lowered:
        return "SEY1/RHD3-like membrane-remodeling GTPase"
    if "cyclin-g-associated kinase" in lowered:
        return "protein kinase"
    if "glycosyltransferase" in lowered:
        return "glycosyltransferase-like enzyme"
    if "mannosyltransferase" in lowered:
        return "mannosyltransferase-like enzyme"
    if "trypsin" in lowered or "serine protease" in lowered or "chymotrypsin" in lowered:
        return "trypsin-like serine protease"
    return description


def infer_function(row: dict[str, object]) -> str:
    pfam = str(row["top_pfam_description"])
    blast = str(row["top_blast_description"])
    signalp_fraction = float(row["signalp_fraction"])
    tm_fraction = float(row["tm_fraction"])
    blast_count = int(row.get("top_blast_count", 0) or 0)
    blast_annotated = int(row.get("blast_annotated_queries", 0) or 0)
    pfam_count = int(row.get("top_pfam_count", 0) or 0)
    pfam_annotated = int(row.get("pfam_annotated_queries", 0) or 0)
    blast_support = blast_count / blast_annotated if blast_annotated else 0.0
    pfam_support = pfam_count / pfam_annotated if pfam_annotated else 0.0

    pfam_label = family_label_from_pfam(pfam, signalp_fraction, tm_fraction)
    blast_label = family_label_from_blast(blast)
    pfam_is_generic = is_generic_pfam_description(pfam)
    blast_is_membrane_compatible = is_membrane_compatible_label(blast_label)
    sparse_consensus_keywords = (
        "sulfotransferase",
        "glycosyltransferase",
        "mannosyltransferase",
        "sulfatase",
        "protease",
    )

    if blast_label and not is_generic_description(blast):
        if blast_count >= 5 and blast_support >= 0.75 and (
            not pfam_label or pfam_is_generic or "domain" in pfam.lower()
        ):
            return f"Likely {blast_label}"
        if tm_fraction >= 0.5 and not blast_is_membrane_compatible and pfam_count < 2:
            blast_label = ""
    if any(keyword in pfam_label.lower() for keyword in sparse_consensus_keywords) and pfam_count >= 1:
        if any(keyword in blast_label.lower() for keyword in sparse_consensus_keywords) and blast_count >= 2:
            return blast_label
        return pfam_label
    if any(keyword in blast_label.lower() for keyword in sparse_consensus_keywords) and blast_count >= 2:
        return blast_label

    if pfam_label and not is_generic_description(pfam):
        if pfam_count >= 3 and (pfam_support >= 0.35 or blast_count < 3 or not blast_label):
            return pfam_label
        if pfam_count >= 2 and pfam_support >= 0.5 and blast_support < 0.5:
            return pfam_label
    if blast_label and not is_generic_description(blast):
        if blast_count >= 3 and blast_support >= 0.45:
            return f"Likely {blast_label}"
        if blast_count >= 2 and blast_support >= 0.6 and not pfam_label and tm_fraction < 0.5:
            return f"Likely {blast_label}"
    if pfam_label and pfam_count >= 2:
        return pfam_label
    if signalp_fraction >= 0.5 and tm_fraction < 0.3:
        return "Secreted protein family with no strong conserved annotation"
    if tm_fraction >= 0.5:
        return "Membrane protein family with no strong conserved annotation"
    if float(row.get("glob_fraction", 0.0) or 0.0) >= 0.8:
        return "Soluble protein family with weak homology support"
    return "Poorly resolved protein family"


def infer_calcification_relevance(row: dict[str, object]) -> tuple[str, str]:
    evidence = " ".join(
        [
            str(row["function_call"]),
            str(row["top_pfam_description"]),
            str(row["top_blast_description"]),
            str(row["top_pfam_summary"]),
            str(row["top_blast_summary"]),
        ]
    ).lower()
    signalp_fraction = float(row["signalp_fraction"])
    tm_fraction = float(row["tm_fraction"])

    if any(keyword in evidence for keyword in DIRECT_CALCIFICATION_KEYWORDS):
        return (
            "likely direct",
            "Contains transporter/ion-homeostasis or carbon-handling signatures that could act directly in calcification chemistry.",
        )
    if signalp_fraction >= 0.5 and tm_fraction < 0.3:
        return (
            "plausible indirect",
            "Secreted or lumenal architecture could fit extracellular matrix, coccolith-associated, or secretory-pathway functions.",
        )
    if tm_fraction >= 0.5:
        return (
            "plausible indirect",
            "Membrane localization makes roles in transport, signaling, trafficking, or vesicle physiology plausible, but the family call is not specific enough for a direct calcification claim.",
        )
    if any(keyword in evidence for keyword in PLAUSIBLE_CALCIFICATION_KEYWORDS):
        return (
            "plausible indirect",
            "Function points to trafficking, cytoskeletal organization, redox control, or secretory biology that can modulate calcification indirectly.",
        )
    if any(keyword in evidence for keyword in LOW_RELEVANCE_KEYWORDS):
        return (
            "low/indirect",
            "Evidence points to general housekeeping rather than a calcification-specific role.",
        )
    return (
        "low/unclear",
        "Available homology and domain evidence does not point to an obvious calcification-related role.",
    )


def build_rows(root: Path) -> list[dict[str, object]]:
    fasta_paths = sorted(root.glob("OG*.fa"))
    fasta_records_by_og: dict[str, list[dict[str, object]]] = {}
    surrogate_to_original: dict[str, str] = {}
    for fasta_path in fasta_paths:
        records = parse_fasta(fasta_path)
        fasta_records_by_og[fasta_path.stem] = records
        for record in records:
            surrogate_to_original[str(record["surrogate_id"])] = str(record["original_id"])

    signalp_by_query: dict[str, dict[str, object]] = {}
    for path in sorted((root / "signalp_results").glob("*.output.json")):
        for row in parse_signalp_json(path):
            original_id = surrogate_to_original.get(str(row["sequence_id"]))
            if original_id is None:
                continue
            signalp_by_query[original_id] = row

    deeptmhmm_by_query: dict[str, dict[str, object]] = {}
    for path in sorted((root / "deeptmhmm_results").glob("*.predicted_topologies.3line")):
        for row in parse_deeptmhmm_3line(path):
            original_id = surrogate_to_original.get(str(row["sequence_id"]))
            if original_id is None:
                continue
            deeptmhmm_by_query[original_id] = row

    rows: list[dict[str, object]] = []
    for fasta_path in fasta_paths:
        orthogroup = fasta_path.stem
        fasta_records = fasta_records_by_og[orthogroup]
        blast_rows = parse_blast_file(root / "blast_out" / f"{orthogroup}.blast.tsv")
        pfam_rows = parse_pfam_tbl(root / "hmm_out" / f"{orthogroup}.pfam.tbl")

        best_blast_by_query: dict[str, dict[str, object]] = {}
        for hit in blast_rows:
            query = str(hit["query"])
            best = best_blast_by_query.get(query)
            if best is None or (float(hit["evalue"]), -float(hit["bitscore"])) < (
                float(best["evalue"]),
                -float(best["bitscore"]),
            ):
                best_blast_by_query[query] = hit

        blast_desc_counter = Counter(
            str(hit["normalized_description"])
            for hit in best_blast_by_query.values()
            if str(hit["normalized_description"])
        )

        pfam_descs_by_query: dict[str, set[str]] = defaultdict(set)
        for domain_row in pfam_rows:
            pfam_descs_by_query[str(domain_row["query"])].add(str(domain_row["description"]))

        pfam_desc_counter = Counter()
        for descriptions in pfam_descs_by_query.values():
            pfam_desc_counter.update(descriptions)

        signalp_positive = 0
        deeptmhmm_labels = Counter()
        for record in fasta_records:
            original_id = str(record["original_id"])
            signalp_row = signalp_by_query.get(original_id)
            if signalp_row and str(signalp_row["prediction"]) != "Other":
                signalp_positive += 1
            deeptmhmm_row = deeptmhmm_by_query.get(original_id)
            if deeptmhmm_row:
                deeptmhmm_labels.update([str(deeptmhmm_row["label"])])

        seq_count = len(fasta_records)
        lengths = [int(record["length"]) for record in fasta_records]
        tm_positive = sum(deeptmhmm_labels.get(label, 0) for label in ("TM", "SP+TM", "BETA"))
        row: dict[str, object] = {
            "orthogroup": orthogroup,
            "sequence_count": seq_count,
            "min_length": min(lengths),
            "median_length": int(statistics.median(lengths)),
            "max_length": max(lengths),
            "blast_annotated_queries": len(best_blast_by_query),
            "pfam_annotated_queries": len(pfam_descs_by_query),
            "signalp_positive": signalp_positive,
            "signalp_fraction": round(signalp_positive / seq_count, 3),
            "deeptmhmm_labels": summarize_counter(deeptmhmm_labels),
            "glob_fraction": round(deeptmhmm_labels.get("GLOB", 0) / seq_count, 3),
            "tm_fraction": round(tm_positive / seq_count, 3),
            "top_blast_description": top_informative_label(blast_desc_counter),
            "top_blast_summary": summarize_counter(blast_desc_counter),
            "median_top_blast_pident": round(median_or_zero(hit["pident"] for hit in best_blast_by_query.values()), 2),
            "median_top_blast_evalue": median_or_zero(hit["evalue"] for hit in best_blast_by_query.values()),
            "top_pfam_description": top_informative_label(pfam_desc_counter),
            "top_pfam_summary": summarize_counter(pfam_desc_counter),
        }
        row["top_blast_count"] = blast_desc_counter.get(str(row["top_blast_description"]), 0)
        row["top_pfam_count"] = pfam_desc_counter.get(str(row["top_pfam_description"]), 0)
        row["function_call"] = infer_function(row)
        relevance, rationale = infer_calcification_relevance(row)
        row["calcification_relevance"] = relevance
        row["calcification_rationale"] = rationale
        rows.append(row)
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "orthogroup",
        "sequence_count",
        "min_length",
        "median_length",
        "max_length",
        "blast_annotated_queries",
        "pfam_annotated_queries",
        "signalp_positive",
        "signalp_fraction",
        "glob_fraction",
        "tm_fraction",
        "deeptmhmm_labels",
        "top_blast_description",
        "top_blast_count",
        "top_blast_summary",
        "median_top_blast_pident",
        "median_top_blast_evalue",
        "top_pfam_description",
        "top_pfam_count",
        "top_pfam_summary",
        "function_call",
        "calcification_relevance",
        "calcification_rationale",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_markdown(path: Path, rows: list[dict[str, object]]) -> None:
    lines = [
        "# Orthogroup Functional Annotation Report",
        "",
        "This report summarizes orthogroup-level evidence from BLAST, PFAM, SignalP, and DeepTMHMM.",
        "",
        "| Orthogroup | Function | Calcification relevance | Key evidence |",
        "| --- | --- | --- | --- |",
    ]
    for row in rows:
        evidence = (
            f"BLAST: {row['top_blast_summary'] or 'n/a'}; "
            f"PFAM: {row['top_pfam_summary'] or 'n/a'}; "
            f"SignalP+: {row['signalp_positive']}/{row['sequence_count']}; "
            f"DeepTMHMM: {row['deeptmhmm_labels'] or 'n/a'}"
        )
        lines.append(
            f"| {row['orthogroup']} | {row['function_call']} | "
            f"{row['calcification_relevance']} | {evidence} |"
        )
        lines.append(
            f"|  |  |  | Rationale: {row['calcification_rationale']} |"
        )
    path.write_text("\n".join(lines) + "\n")


def main(argv: list[str]) -> int:
    root = Path(argv[1]).resolve() if len(argv) > 1 else Path(__file__).resolve().parent
    rows = build_rows(root)
    write_tsv(root / "orthogroup_annotation_evidence.tsv", rows)
    write_markdown(root / "orthogroup_annotation_report.md", rows)
    print(f"Wrote {len(rows)} orthogroup summaries to {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
