#!/usr/bin/env python3
"""Summarize BLAST, Pfam (HMM), DeepTMHMM, SignalP evidence per orthogroup."""
from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BLAST_DIR = ROOT / "blast_out"
HMM_DIR = ROOT / "hmm_out"
TM_DIR = ROOT / "deeptmhmm_results"
SP_DIR = ROOT / "signalp_results"

OG_RE = re.compile(r"OG\d{7}")


def list_ogs() -> list[str]:
    return sorted({p.stem for p in ROOT.glob("OG*.fa")})


def count_fa(og: str) -> int:
    fa = ROOT / f"{og}.fa"
    n = 0
    with fa.open() as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def parse_blast(og: str) -> dict:
    """Return best-hit stats and top hit descriptions for an OG."""
    path = BLAST_DIR / f"{og}.blast.tsv"
    if not path.exists() or path.stat().st_size == 0:
        return {"n_queries_with_hit": 0, "top_subjects": [], "best_descs": []}

    # Per query, keep best (lowest evalue, highest bitscore).
    best_per_query: dict[str, tuple[float, float, str, str]] = {}
    with path.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            q, subj, desc, pid, length, evalue, bitscore = (
                parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]
            )
            try:
                ev = float(evalue)
                bs = float(bitscore)
            except ValueError:
                continue
            cur = best_per_query.get(q)
            if cur is None or (ev, -bs) < (cur[0], -cur[1]):
                best_per_query[q] = (ev, bs, subj, desc)

    subj_counter: Counter[str] = Counter()
    desc_counter: Counter[str] = Counter()
    for ev, bs, subj, desc in best_per_query.values():
        subj_counter[subj] += 1
        # Strip the leading redundant "sp|..|.. " prefix if present to get description phrase only
        clean = desc
        m = re.match(r"^sp\|[^|]+\|\S+\s+(.*?)\s+OS=", desc)
        if m:
            clean = m.group(1)
        desc_counter[clean] += 1

    return {
        "n_queries_with_hit": len(best_per_query),
        "top_subjects": subj_counter.most_common(5),
        "best_descs": desc_counter.most_common(5),
    }


def parse_pfam(og: str) -> dict:
    """Count Pfam domain hits for proteins in the OG."""
    path = HMM_DIR / f"{og}.pfam.tbl"
    dom_counter: Counter[str] = Counter()
    dom_desc: dict[str, str] = {}
    n_queries_with_domain: set[str] = set()
    if not path.exists():
        return {
            "n_queries_with_domain": 0,
            "top_domains": [],
            "dom_desc": {},
        }
    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            # hmmscan tblout: target_name accession query_name ...
            toks = line.split(None, 18)
            if len(toks) < 19:
                continue
            tname, tacc, qname = toks[0], toks[1], toks[2]
            # E-value full is toks[4]; filter on e-value < 1e-3
            try:
                evalue = float(toks[4])
            except ValueError:
                continue
            if evalue > 1e-3:
                continue
            desc = toks[18].strip()
            dom_counter[tname] += 1
            dom_desc[tname] = desc
            n_queries_with_domain.add(qname)
    return {
        "n_queries_with_domain": len(n_queries_with_domain),
        "top_domains": dom_counter.most_common(8),
        "dom_desc": dom_desc,
    }


def parse_deeptmhmm() -> dict[str, dict]:
    """Return per-OG summary of topology types and TM counts."""
    per_protein: dict[str, str] = {}
    for f in TM_DIR.glob("*.predicted_topologies.3line"):
        with f.open() as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                _seq = fh.readline()
                _top = fh.readline()
                if not header.startswith(">"):
                    continue
                # >OG0000049__000001 | GLOB
                m = re.match(r">(OG\d+__\d+)\s*\|\s*(\S+)", header.strip())
                if not m:
                    continue
                per_protein[m.group(1)] = m.group(2)

    # Summarize
    per_og: dict[str, Counter] = defaultdict(Counter)
    for pid, ttype in per_protein.items():
        og = pid.split("__")[0]
        per_og[og][ttype] += 1
    return {og: dict(c) for og, c in per_og.items()}


def parse_signalp() -> dict[str, dict]:
    per_protein: dict[str, str] = {}
    for f in SP_DIR.glob("signalp_batch_*.output.json"):
        data = json.loads(f.read_text())
        for pid, info in data.get("SEQUENCES", {}).items():
            per_protein[pid] = info.get("Prediction", "Other")
    per_og: dict[str, Counter] = defaultdict(Counter)
    for pid, pred in per_protein.items():
        og = pid.split("__")[0]
        per_og[og][pred] += 1
    return {og: dict(c) for og, c in per_og.items()}


def main() -> None:
    ogs = list_ogs()
    tm_summary = parse_deeptmhmm()
    sp_summary = parse_signalp()

    out_path = ROOT / "annotation_summary.tsv"
    cols = [
        "OG",
        "n_sequences",
        "blast_n_with_hit",
        "top_blast_descriptions",
        "pfam_n_with_domain",
        "top_pfam_domains",
        "signalp_counts",
        "deeptmhmm_counts",
    ]
    with out_path.open("w") as out:
        out.write("\t".join(cols) + "\n")
        for og in ogs:
            n_seq = count_fa(og)
            b = parse_blast(og)
            p = parse_pfam(og)
            sp = sp_summary.get(og, {})
            tm = tm_summary.get(og, {})

            blast_desc_str = "; ".join(f"{d}({n})" for d, n in b["best_descs"])
            pfam_str = "; ".join(
                f"{name}({n})[{p['dom_desc'].get(name,'')}]"
                for name, n in p["top_domains"]
            )
            sp_str = ", ".join(f"{k}:{v}" for k, v in sp.items())
            tm_str = ", ".join(f"{k}:{v}" for k, v in tm.items())

            row = [
                og,
                str(n_seq),
                str(b["n_queries_with_hit"]),
                blast_desc_str,
                str(p["n_queries_with_domain"]),
                pfam_str,
                sp_str,
                tm_str,
            ]
            out.write("\t".join(row) + "\n")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
