from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


NOISY_BLAST_TERMS = {
    "hypothetical protein",
    "uncharacterized protein",
    "predicted protein",
    "unknown protein",
    "fragment",
}

DIRECT_CALCIFICATION_KEYWORDS = {
    "carbonic anhydrase": 4,
    "bicarbonate": 4,
    "calcium": 4,
    "ca2+": 4,
    "shell matrix": 4,
    "adhesive plaque": 4,
    "extracellular": 1,
    "laminin": 2,
    "collagen": 2,
    "cadherin": 2,
    "von willebrand": 2,
    "egf": 1,
    "mucin": 2,
    "sushi": 2,
    "lectin": 2,
    "tyrosinase": 3,
    "peroxidase": 3,
    "multicopper": 3,
    "chitin": 3,
    "glycosyltransferase": 3,
    "mannosyltransferase": 2,
    "exostosin": 3,
    "sulfatase": 2,
    "sulfotransferase": 3,
    "trypsin": 2,
    "polysaccharide": 2,
    "anion transporter": 3,
    "solute carrier": 2,
    "channel": 2,
    "transporter": 2,
    "v-type proton": 2,
    "atpase": 1,
    "secreted": 1,
    "acidic": 1,
    "cysteine-rich": 1,
}

INDIRECT_CALCIFICATION_KEYWORDS = {
    "vesicle": 3,
    "golgi": 3,
    "snare": 3,
    "trafficking": 3,
    "exocytosis": 3,
    "cytoskeleton": 2,
    "actin": 2,
    "microtubule": 2,
    "coiled-coil": 1,
    "cilia": 2,
    "adhesion": 2,
    "integrin": 2,
    "kinase": 1,
    "g protein": 1,
    "receptor": 1,
    "membrane": 1,
    "cell wall": 2,
}

UNLIKELY_CALCIFICATION_KEYWORDS = {
    "rna polymerase": 4,
    "transcription": 3,
    "ribosomal": 4,
    "ribosome": 4,
    "replication": 3,
    "repair": 3,
    "endonuclease": 3,
    "helicase": 3,
    "proteasome": 3,
    "splicing": 3,
    "translation": 4,
    "mitotic": 2,
    "checkpoint": 2,
    "nuclear": 2,
    "dna-directed": 4,
}

FUNCTION_RULES = [
    (("carbonic anhydrase",), "carbonic anhydrase"),
    (("bicarbonate transporter", "anion transporter"), "anion/bicarbonate transporter"),
    (("calcium channel", "calcium-transporting", "calcium transporter"), "calcium transporter/channel"),
    (("fk506-binding protein 15", "fkbp15"), "FKBP15-like coiled-coil scaffold protein"),
    (("coiled-coil domain-containing protein 39", "ccdc39", "odf2"), "coiled-coil scaffold or cytoskeletal protein"),
    (("warthog protein", "hedgehog", "hint module"), "hedgehog/warthog-like signaling protein"),
    (("rcc1-like",), "RCC1-like repeat protein"),
    (("protein kinase", "serine/threonine kinase", "cyclin-g-associated kinase"), "protein kinase"),
    (("tyrosinase", "multicopper oxidase"), "oxidase involved in extracellular matrix hardening"),
    (("peroxidase",), "peroxidase"),
    (("chitin synthase",), "chitin synthase"),
    (("chitinase",), "chitin remodeling enzyme"),
    (("sulfatase", "arylsulfatase"), "arylsulfatase / extracellular sulfatase"),
    (("sulfotransferase",), "glycan sulfotransferase"),
    (("fg-gap", "integrin alpha"), "integrin alpha-like cell-surface receptor"),
    (("adhesive plaque", "collagen", "laminin", "von willebrand", "cadherin", "mucin", "sushi", "fasciclin", "egf"), "extracellular matrix or adhesion protein"),
    (("glycosyl transferase family 8", "muci70", "tod1"), "matrix polysaccharide glycosyltransferase"),
    (("mannosyltransferase",), "alpha-1,3-mannosyltransferase-like glycosyltransferase"),
    (("exostosin",), "exostosin-like transmembrane glycosyltransferase"),
    (("glycosyltransferase", "polysaccharide"), "glycosyltransferase involved in extracellular polysaccharide synthesis"),
    (("trypsin", "serine protease"), "secreted serine protease"),
    (("lectin",), "lectin-like extracellular protein"),
    (("prima1", "shoct"), "short membrane-anchored protein"),
    (("sur7/pali", "sur7"), "SUR7/PalI-like membrane organizer"),
    (("rhd3", "sey1"), "ER membrane fusion dynamin-like GTPase"),
    (("fad binding domain", "monooxygenase", "oxidoreductase fad-binding"), "FAD-dependent oxidoreductase"),
    (("amidase", "amidotransferase"), "amidase or amidotransferase-like enzyme"),
    (("rna polymerase ii subunit rpb1", "polr2a"), "RNA polymerase II largest subunit"),
    (("endonuclease iii", "nthl1"), "DNA repair glycosylase"),
]


def parse_float(text: str) -> float:
    return float(text)


def safe_mean(values: Iterable[float]) -> float | None:
    values = list(values)
    if not values:
        return None
    return sum(values) / len(values)


def safe_median(values: Iterable[float]) -> float | None:
    values = list(values)
    if not values:
        return None
    return statistics.median(values)


def format_number(value: float | None, digits: int = 3) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def format_evalue(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.1e}"


def canonical_sequence_id(og_id: str, index: int) -> str:
    return f"{og_id}__{index:06d}"


def aa_fraction(sequence: str, residues: str) -> float:
    if not sequence:
        return 0.0
    wanted = sum(sequence.count(residue) for residue in residues)
    return wanted / len(sequence)


def low_complexity_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    counts = Counter(sequence)
    return counts.most_common(1)[0][1] / len(sequence)


def parse_fasta(path: Path, og_id: str) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    header: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal header, seq_parts
        if header is None:
            return
        original_id = header.split()[0]
        index = len(records) + 1
        sequence = "".join(seq_parts).replace("*", "")
        records.append(
            {
                "generated_id": canonical_sequence_id(og_id, index),
                "original_id": original_id,
                "original_header": header,
                "length": len(sequence),
                "acidic_fraction": aa_fraction(sequence, "DE"),
                "basic_fraction": aa_fraction(sequence, "KRH"),
                "cysteine_fraction": aa_fraction(sequence, "C"),
                "glycine_fraction": aa_fraction(sequence, "G"),
                "low_complexity_fraction": low_complexity_fraction(sequence),
            }
        )
        header = None
        seq_parts = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:]
            else:
                seq_parts.append(line)
    flush()
    return records


def parse_signalp_cs_pos(text: str) -> int | None:
    if not text:
        return None
    match = re.search(r"between pos\.\s+(\d+)\s+and\s+(\d+)", text)
    if not match:
        return None
    return int(match.group(1))


def load_signalp_predictions(signalp_dir: Path) -> dict[str, dict[str, object]]:
    predictions: dict[str, dict[str, object]] = {}
    for json_path in sorted(signalp_dir.glob("signalp_batch_*.output.json")):
        with json_path.open() as handle:
            payload = json.load(handle)
        for seq_id, record in payload.get("SEQUENCES", {}).items():
            protein_types = record.get("Protein_types", [])
            likelihoods = record.get("Likelihood", [])
            score_by_type = {
                protein_type: float(score)
                for protein_type, score in zip(protein_types, likelihoods)
            }
            predictions[seq_id] = {
                "prediction": record.get("Prediction", "Other"),
                "sp_probability": score_by_type.get("Signal Peptide (Sec/SPI)", 0.0),
                "cleavage_site_end": parse_signalp_cs_pos(record.get("CS_pos", "")),
            }
    return predictions


def count_runs(text: str, residue: str) -> int:
    count = 0
    in_run = False
    for char in text:
        if char == residue:
            if not in_run:
                count += 1
                in_run = True
        else:
            in_run = False
    return count


def load_deeptmhmm_predictions(deeptmhmm_dir: Path) -> dict[str, dict[str, object]]:
    predictions: dict[str, dict[str, object]] = {}
    for path in sorted(deeptmhmm_dir.glob("*.predicted_topologies.3line")):
        with path.open() as handle:
            lines = [line.rstrip("\n") for line in handle if line.strip()]
        for index in range(0, len(lines), 3):
            if index + 2 >= len(lines):
                continue
            header, _, topology = lines[index : index + 3]
            match = re.match(r"^>(\S+)\s+\|\s+(\S+)", header)
            if not match:
                continue
            seq_id, predicted_class = match.groups()
            predictions[seq_id] = {
                "class": predicted_class,
                "tm_count": count_runs(topology, "M"),
            }
    return predictions


def parse_pfam_tbl(path: Path) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split(maxsplit=18)
            if len(fields) < 19:
                continue
            hits.append(
                {
                    "domain_name": fields[0],
                    "domain_accession": fields[1],
                    "query_name": fields[2],
                    "full_evalue": parse_float(fields[4]),
                    "domain_evalue": parse_float(fields[7]),
                    "description": fields[18],
                }
            )
    return hits


def clean_blast_description(description: str) -> str:
    cleaned = description.split(" OS=")[0].strip()
    cleaned = re.sub(r"^(?:sp|tr|ref)\|[^|]+\|\S+\s+", "", cleaned)
    cleaned = re.sub(r"\s+", " ", cleaned)
    return cleaned


def parse_blast_tsv(path: Path) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or len(row) < 7:
                continue
            hits.append(
                {
                    "query_name": row[0],
                    "subject_id": row[1],
                    "subject_description": clean_blast_description(row[2]),
                    "percent_identity": float(row[3]),
                    "alignment_length": int(float(row[4])),
                    "evalue": float(row[5]),
                    "bitscore": float(row[6]),
                }
            )
    hits.sort(key=lambda hit: (hit["query_name"], hit["evalue"], -hit["bitscore"]))
    return hits


def filter_pfam_hits(hits: list[dict[str, object]]) -> list[dict[str, object]]:
    filtered = []
    for hit in hits:
        if float(hit["full_evalue"]) > 1e-4:
            continue
        if float(hit["domain_evalue"]) > 1e-3:
            continue
        filtered.append(hit)
    return filtered


def filter_reliable_blast_hits(hits: list[dict[str, object]]) -> list[dict[str, object]]:
    reliable: list[dict[str, object]] = []
    best_per_query: dict[str, dict[str, object]] = {}
    for hit in hits:
        query_name = str(hit["query_name"])
        best_per_query.setdefault(query_name, hit)
    for hit in best_per_query.values():
        evalue = float(hit["evalue"])
        alignment_length = int(hit["alignment_length"])
        if evalue <= 1e-10:
            reliable.append(hit)
        elif evalue <= 1e-4 and alignment_length >= 80:
            reliable.append(hit)
    reliable.sort(key=lambda hit: (hit["query_name"], hit["evalue"], -hit["bitscore"]))
    return reliable


def summarize_pfam_hits(hits: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], dict[str, object]] = {}
    seen_pairs: set[tuple[str, str]] = set()
    for hit in hits:
        key = (str(hit["domain_name"]), str(hit["description"]))
        query_name = str(hit["query_name"])
        entry = grouped.setdefault(
            key,
            {
                "domain_name": hit["domain_name"],
                "description": hit["description"],
                "best_evalue": hit["domain_evalue"],
                "queries": set(),
            },
        )
        entry["best_evalue"] = min(float(entry["best_evalue"]), float(hit["domain_evalue"]))
        pair_key = (str(hit["domain_name"]), query_name)
        if pair_key not in seen_pairs:
            entry["queries"].add(query_name)
            seen_pairs.add(pair_key)
    summaries = []
    for entry in grouped.values():
        summaries.append(
            {
                "domain_name": entry["domain_name"],
                "description": entry["description"],
                "count": len(entry["queries"]),
                "best_evalue": entry["best_evalue"],
            }
        )
    summaries.sort(key=lambda item: (-int(item["count"]), float(item["best_evalue"])))
    return summaries


def summarize_blast_hits(hits: list[dict[str, object]]) -> list[dict[str, object]]:
    best_per_query: dict[str, dict[str, object]] = {}
    for hit in hits:
        query_name = str(hit["query_name"])
        best_per_query.setdefault(query_name, hit)
    grouped: dict[str, dict[str, object]] = {}
    for hit in best_per_query.values():
        description = str(hit["subject_description"])
        entry = grouped.setdefault(
            description,
            {"description": description, "count": 0, "best_evalue": hit["evalue"]},
        )
        entry["count"] += 1
        entry["best_evalue"] = min(float(entry["best_evalue"]), float(hit["evalue"]))
    summaries = list(grouped.values())
    summaries.sort(
        key=lambda item: (
            description_is_noisy(str(item["description"])),
            -int(item["count"]),
            float(item["best_evalue"]),
        )
    )
    return summaries


def description_is_noisy(description: str) -> bool:
    normalized = description.lower()
    return any(term in normalized for term in NOISY_BLAST_TERMS)


def format_summary_items(items: list[dict[str, object]], label_key: str, limit: int = 4) -> str:
    pieces = []
    for item in items[:limit]:
        pieces.append(
            f"{item[label_key]} (n={item['count']}, bestE={format_evalue(float(item['best_evalue']))})"
        )
    return "; ".join(pieces)


def infer_function_label(
    pfam_summary: list[dict[str, object]],
    blast_summary: list[dict[str, object]],
    signal_fraction: float,
    secreted_acidic_mean: float | None,
    secreted_cysteine_mean: float | None,
) -> str:
    evidence_parts: list[str] = []
    for item in pfam_summary[:4]:
        evidence_parts.append(str(item["domain_name"]).lower())
        evidence_parts.append(str(item["description"]).lower())
    for item in blast_summary[:4]:
        evidence_parts.append(str(item["description"]).lower())
    evidence_text = " ".join(evidence_parts)

    for keywords, label in FUNCTION_RULES:
        if any(keyword in evidence_text for keyword in keywords):
            return label

    if signal_fraction >= 0.5 and (secreted_acidic_mean or 0.0) >= 0.13:
        return "secreted acidic protein"
    if signal_fraction >= 0.5 and (secreted_cysteine_mean or 0.0) >= 0.04:
        return "secreted cysteine-rich protein"
    if pfam_summary:
        return str(pfam_summary[0]["description"])
    if blast_summary:
        return str(blast_summary[0]["description"])
    if signal_fraction >= 0.5:
        return "secreted protein with no strong domain assignment"
    return "uncharacterized protein family"


def score_keyword_set(text: str, weights: dict[str, int]) -> int:
    score = 0
    for keyword, weight in weights.items():
        if keyword in text:
            score += weight
    return score


def assign_calcification_relevance(
    function_label: str,
    pfam_summary: list[dict[str, object]],
    blast_summary: list[dict[str, object]],
    seq_count: int,
    signal_count: int,
    tm_count: int,
    beta_count: int,
    sp_tm_count: int,
    secreted_acidic_mean: float | None,
    secreted_cysteine_mean: float | None,
) -> tuple[str, str]:
    evidence_parts = [function_label.lower()]
    for item in pfam_summary[:4]:
        evidence_parts.append(str(item["domain_name"]).lower())
        evidence_parts.append(str(item["description"]).lower())
    for item in blast_summary[:4]:
        evidence_parts.append(str(item["description"]).lower())
    evidence_text = " ".join(evidence_parts)
    strong_direct_feature = any(
        keyword in evidence_text
        for keyword in (
            "carbonic anhydrase",
            "bicarbonate",
            "calcium",
            "adhesive plaque",
            "collagen",
            "laminin",
            "glycosyltransferase",
            "muci70",
            "tod1",
            "exostosin",
            "mannosyltransferase",
            "sulfatase",
            "sulfotransferase",
            "trypsin",
            "chitin",
            "tyrosinase",
            "peroxidase",
        )
    )

    direct_score = score_keyword_set(evidence_text, DIRECT_CALCIFICATION_KEYWORDS)
    indirect_score = score_keyword_set(evidence_text, INDIRECT_CALCIFICATION_KEYWORDS)
    unlikely_score = score_keyword_set(evidence_text, UNLIKELY_CALCIFICATION_KEYWORDS)

    signal_fraction = signal_count / seq_count if seq_count else 0.0
    tm_fraction = tm_count / seq_count if seq_count else 0.0
    beta_fraction = beta_count / seq_count if seq_count else 0.0
    sp_tm_fraction = sp_tm_count / seq_count if seq_count else 0.0
    secreted_non_tm_fraction = (signal_count - sp_tm_count) / seq_count if seq_count else 0.0

    if signal_fraction >= 0.5 and tm_fraction <= 0.2:
        direct_score += 2
    if beta_fraction >= 0.3 or sp_tm_fraction >= 0.3:
        indirect_score += 1
    if tm_fraction >= 0.5 and not any(
        keyword in evidence_text for keyword in ("transporter", "channel", "pump")
    ):
        direct_score = max(0, direct_score - 1)
        indirect_score += 1
    if "kinase" in evidence_text and signal_fraction == 0 and tm_fraction == 0:
        direct_score = max(0, direct_score - 2)
        indirect_score += 1
    if (secreted_acidic_mean or 0.0) >= 0.13:
        direct_score += 2
    if (secreted_cysteine_mean or 0.0) >= 0.04:
        direct_score += 1
    if signal_fraction == 0 and tm_fraction == 0:
        unlikely_score += 1

    reasons = []
    if signal_count:
        reasons.append(f"{signal_count}/{seq_count} sequences have a signal peptide")
    if tm_count:
        reasons.append(f"{tm_count}/{seq_count} sequences are predicted transmembrane")
    if beta_count:
        reasons.append(f"{beta_count}/{seq_count} sequences are predicted beta-barrel")
    if pfam_summary:
        reasons.append(f"PFAM consensus: {pfam_summary[0]['description']}")
    elif blast_summary:
        reasons.append(f"BLAST consensus: {blast_summary[0]['description']}")
    if (secreted_acidic_mean or 0.0) >= 0.13:
        reasons.append("secreted members are acidic")
    if (secreted_cysteine_mean or 0.0) >= 0.04:
        reasons.append("secreted members are cysteine-rich")

    if (
        direct_score >= 6
        and strong_direct_feature
        and unlikely_score <= 2
        and (secreted_non_tm_fraction >= 0.2 or any(
            keyword in evidence_text for keyword in ("transporter", "channel", "pump", "glycosyltransferase")
        ))
    ):
        relevance = "high"
        reasons.append("features fit a plausible direct role in biomineral deposition or ion handling")
    elif direct_score >= 3 or indirect_score >= 4:
        relevance = "moderate"
        reasons.append("features fit extracellular, transport, or trafficking roles that could influence calcification indirectly")
    elif unlikely_score >= 4:
        relevance = "unlikely"
        reasons.append("evidence points to a housekeeping intracellular role with no obvious calcification link")
    else:
        relevance = "low"
        reasons.append("current evidence does not suggest a strong calcification-specific function")
    return relevance, "; ".join(reasons)


def summarize_orthogroup(
    root: Path,
    og_id: str,
    signalp_predictions: dict[str, dict[str, object]],
    deeptmhmm_predictions: dict[str, dict[str, object]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    fasta_path = root / f"{og_id}.fa"
    blast_path = root / "blast_out" / f"{og_id}.blast.tsv"
    pfam_path = root / "hmm_out" / f"{og_id}.pfam.tbl"

    records = parse_fasta(fasta_path, og_id)
    by_original_id = {str(record["original_id"]): record for record in records}

    blast_hits = parse_blast_tsv(blast_path)
    pfam_hits = parse_pfam_tbl(pfam_path)
    filtered_pfam_hits = filter_pfam_hits(pfam_hits)
    filtered_blast_hits = filter_reliable_blast_hits(blast_hits)

    best_blast_by_query: dict[str, dict[str, object]] = {}
    for hit in filtered_blast_hits:
        query_name = str(hit["query_name"])
        best_blast_by_query.setdefault(query_name, hit)

    pfam_by_query: dict[str, list[dict[str, object]]] = defaultdict(list)
    for hit in filtered_pfam_hits:
        pfam_by_query[str(hit["query_name"])].append(hit)

    for record in records:
        generated_id = str(record["generated_id"])
        original_id = str(record["original_id"])
        signalp_record = signalp_predictions.get(generated_id, {})
        tm_record = deeptmhmm_predictions.get(generated_id, {})
        blast_record = best_blast_by_query.get(original_id)
        record["signalp_prediction"] = signalp_record.get("prediction", "Other")
        record["signalp_probability"] = signalp_record.get("sp_probability", 0.0)
        record["cleavage_site_end"] = signalp_record.get("cleavage_site_end")
        record["deeptmhmm_class"] = tm_record.get("class", "")
        record["tm_count"] = tm_record.get("tm_count", 0)
        record["blast_top_description"] = (
            blast_record["subject_description"] if blast_record else ""
        )
        record["blast_top_evalue"] = blast_record["evalue"] if blast_record else None
        record["pfam_domains"] = "; ".join(
            sorted({str(hit["domain_name"]) for hit in pfam_by_query.get(original_id, [])})
        )

    pfam_summary = summarize_pfam_hits(filtered_pfam_hits)
    blast_summary = summarize_blast_hits(filtered_blast_hits)

    lengths = [int(record["length"]) for record in records]
    signal_count = sum(record["signalp_prediction"] != "Other" for record in records)
    tm_positive_count = sum(int(record["tm_count"]) > 0 for record in records)
    beta_count = sum(record["deeptmhmm_class"] == "BETA" for record in records)
    sp_tm_count = sum(record["deeptmhmm_class"] == "SP+TM" for record in records)
    secreted_non_tm = [
        record
        for record in records
        if record["signalp_prediction"] != "Other" and int(record["tm_count"]) == 0
    ]
    secreted_acidic_mean = safe_mean(
        float(record["acidic_fraction"]) for record in secreted_non_tm
    )
    secreted_cysteine_mean = safe_mean(
        float(record["cysteine_fraction"]) for record in secreted_non_tm
    )

    function_label = infer_function_label(
        pfam_summary,
        blast_summary,
        signal_count / len(records) if records else 0.0,
        secreted_acidic_mean,
        secreted_cysteine_mean,
    )
    calcification_relevance, rationale = assign_calcification_relevance(
        function_label,
        pfam_summary,
        blast_summary,
        len(records),
        signal_count,
        tm_positive_count,
        beta_count,
        sp_tm_count,
        secreted_acidic_mean,
        secreted_cysteine_mean,
    )

    summary = {
        "orthogroup": og_id,
        "seq_count": len(records),
        "median_length": int(safe_median(lengths) or 0),
        "mean_length": format_number(safe_mean(lengths), digits=1),
        "signal_peptide_count": signal_count,
        "signal_peptide_fraction": format_number(signal_count / len(records) if records else 0.0),
        "tm_protein_count": tm_positive_count,
        "beta_barrel_count": beta_count,
        "sp_tm_count": sp_tm_count,
        "secreted_non_tm_count": len(secreted_non_tm),
        "secreted_acidic_fraction_mean": format_number(secreted_acidic_mean),
        "secreted_cysteine_fraction_mean": format_number(secreted_cysteine_mean),
        "top_pfam_domains": format_summary_items(pfam_summary, "domain_name"),
        "top_pfam_descriptions": format_summary_items(pfam_summary, "description"),
        "top_blast_consensus": format_summary_items(blast_summary, "description"),
        "inferred_function": function_label,
        "calcification_relevance": calcification_relevance,
        "calcification_rationale": rationale,
    }
    return summary, records


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown_report(path: Path, summary_rows: list[dict[str, object]]) -> None:
    high = [row for row in summary_rows if row["calcification_relevance"] == "high"]
    moderate = [row for row in summary_rows if row["calcification_relevance"] == "moderate"]
    unlikely = [row for row in summary_rows if row["calcification_relevance"] == "unlikely"]

    with path.open("w") as handle:
        handle.write("# Orthogroup Function And Calcification Relevance Report\n\n")
        handle.write("This report combines BLAST, PFAM, SignalP, and DeepTMHMM evidence.\n")
        handle.write("Additional features generated here include sequence length statistics plus acidic and cysteine composition for secreted non-membrane members.\n\n")

        def write_section(title: str, rows: list[dict[str, object]]) -> None:
            handle.write(f"## {title}\n\n")
            if not rows:
                handle.write("None.\n\n")
                return
            for row in rows:
                handle.write(f"### {row['orthogroup']}\n")
                handle.write(f"- Function: {row['inferred_function']}\n")
                handle.write(f"- Relevance: {row['calcification_relevance']}\n")
                handle.write(f"- PFAM: {row['top_pfam_descriptions'] or 'none'}\n")
                handle.write(f"- BLAST: {row['top_blast_consensus'] or 'none'}\n")
                handle.write(f"- Rationale: {row['calcification_rationale']}\n\n")

        write_section("Highest-Priority Candidates", high)
        write_section("Moderate Candidates", moderate)
        write_section("Likely Low Or No Direct Relevance", unlikely)

        low_rows = [
            row
            for row in summary_rows
            if row["calcification_relevance"] not in {"high", "moderate", "unlikely"}
        ]
        write_section("Low-Confidence Or Ambiguous Families", low_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--outdir", type=Path, default=None)
    args = parser.parse_args()

    root = args.root.resolve()
    outdir = args.outdir.resolve() if args.outdir else root / "annotation_results"
    outdir.mkdir(parents=True, exist_ok=True)

    signalp_predictions = load_signalp_predictions(root / "signalp_results")
    deeptmhmm_predictions = load_deeptmhmm_predictions(root / "deeptmhmm_results")

    orthogroups = sorted(path.stem for path in root.glob("OG*.fa"))
    summary_rows: list[dict[str, object]] = []
    sequence_rows: list[dict[str, object]] = []

    for og_id in orthogroups:
        summary, records = summarize_orthogroup(
            root,
            og_id,
            signalp_predictions,
            deeptmhmm_predictions,
        )
        summary_rows.append(summary)
        for record in records:
            sequence_rows.append(
                {
                    "orthogroup": og_id,
                    "generated_id": record["generated_id"],
                    "original_id": record["original_id"],
                    "length": record["length"],
                    "acidic_fraction": format_number(float(record["acidic_fraction"])),
                    "cysteine_fraction": format_number(float(record["cysteine_fraction"])),
                    "low_complexity_fraction": format_number(float(record["low_complexity_fraction"])),
                    "signalp_prediction": record["signalp_prediction"],
                    "signalp_probability": format_number(float(record["signalp_probability"])),
                    "deeptmhmm_class": record["deeptmhmm_class"],
                    "tm_count": record["tm_count"],
                    "pfam_domains": record["pfam_domains"],
                    "blast_top_description": record["blast_top_description"],
                    "blast_top_evalue": format_evalue(record["blast_top_evalue"]),
                }
            )

    write_tsv(outdir / "orthogroup_annotation_summary.tsv", summary_rows)
    write_tsv(outdir / "orthogroup_sequence_evidence.tsv", sequence_rows)
    write_markdown_report(outdir / "orthogroup_annotation_report.md", summary_rows)


if __name__ == "__main__":
    main()
