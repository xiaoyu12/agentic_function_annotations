#!/usr/bin/env python3
"""Aggregate annotation evidence per orthogroup.

Produces, per OG:
  - total sequence count
  - species represented
  - top BLAST hit (by bitscore) + its SwissProt entry + % of sequences with any hit
  - top PFAM domains (by occurrence) with E-value range
  - # sequences with a signal peptide (SignalP-6.0)
  - # sequences with at least one predicted TM region (DeepTMHMM)
"""
from __future__ import annotations
import glob
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path("/Users/xiaoyu/workspace/Orthogroups.calcifying_loose_fastas")


def parse_fa(path: Path):
    headers = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


def parse_blast(path: Path):
    """Return list of dicts {query, subject, desc, pident, aln, evalue, bitscore}."""
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            q, s, desc, pid, aln, ev, bs = parts[:7]
            try:
                bs = float(bs)
                ev = float(ev)
                pid = float(pid)
                aln = int(aln)
            except ValueError:
                continue
            rows.append({
                "q": q, "s": s, "desc": desc,
                "pident": pid, "aln": aln,
                "evalue": ev, "bitscore": bs,
            })
    return rows


def parse_pfam(path: Path):
    """Parse hmmsearch --tblout.  Return list of (target_name, acc, query, evalue, score, desc)."""
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split(None, 18)
            if len(parts) < 19:
                continue
            target = parts[0]
            acc = parts[1]
            query = parts[2]
            try:
                evalue = float(parts[4])
                score = float(parts[5])
            except ValueError:
                continue
            desc = parts[18].rstrip("\n")
            rows.append({
                "target": target, "acc": acc,
                "query": query, "evalue": evalue,
                "score": score, "desc": desc,
            })
    return rows


def collect_signalp():
    """Map seq_id -> cleavage site score (presence)."""
    sig = set()
    sig_details = {}  # id -> (start, end, score)
    for gff in sorted(ROOT.glob("signalp_results/*/output.gff3")):
        with open(gff) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6 or parts[2] != "signal_peptide":
                    continue
                sid = parts[0]
                start, end = int(parts[3]), int(parts[4])
                score = float(parts[5])
                sig.add(sid)
                sig_details[sid] = (start, end, score)
    return sig, sig_details


def collect_tm():
    """Map seq_id -> # TM regions."""
    tm_counts = {}
    for fp in sorted(ROOT.glob("deeptmhmm_results/*.TMRs.gff3")):
        with open(fp) as f:
            cur_id = None
            cur_len = None
            cur_tmrs = None
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("# ") and "Length" in line:
                    # "# OG0000049__000001 Length: 901"
                    m = re.match(r"#\s+(\S+)\s+Length:\s+(\d+)", line)
                    if m:
                        cur_id = m.group(1)
                        cur_len = int(m.group(2))
                elif line.startswith("# ") and "Number of predicted TMRs" in line:
                    m = re.match(r"#\s+(\S+)\s+Number of predicted TMRs:\s+(\d+)", line)
                    if m:
                        cur_id = m.group(1)
                        cur_tmrs = int(m.group(2))
                        tm_counts[cur_id] = cur_tmrs
    return tm_counts


def og_of(seq_id: str):
    """Map any seq_id back to its OG id.

    OG___NNNN sequences are obviously OG-prefixed. For JGI headers in blast/hmm,
    we use the inverse map built from fasta files.
    """
    m = re.match(r"(OG\d+)__", seq_id)
    if m:
        return m.group(1)
    return None


def build_header_to_og():
    mapping = {}
    for fa in sorted(ROOT.glob("OG*.fa")):
        og = fa.stem
        for h in parse_fa(fa):
            mapping[h] = og
    return mapping


def main():
    hdr2og = build_header_to_og()
    ogs = sorted({og for og in hdr2og.values()})
    og_seqs = defaultdict(list)
    for h, og in hdr2og.items():
        og_seqs[og].append(h)

    # Species per OG (from JGI headers)
    og_species = defaultdict(Counter)
    for h, og in hdr2og.items():
        sp = h.split("|")[1] if h.startswith("jgi|") else "?"
        og_species[og][sp] += 1

    # BLAST per OG
    og_blast = defaultdict(list)
    og_blast_queries = defaultdict(set)
    for og in ogs:
        rows = parse_blast(ROOT / "blast_out" / f"{og}.blast.tsv")
        for r in rows:
            og_blast[og].append(r)
            og_blast_queries[og].add(r["q"])

    # PFAM per OG
    og_pfam = defaultdict(list)
    og_pfam_queries = defaultdict(set)
    for og in ogs:
        rows = parse_pfam(ROOT / "hmm_out" / f"{og}.pfam.tbl")
        for r in rows:
            og_pfam[og].append(r)
            og_pfam_queries[og].add(r["query"])

    # SignalP + TM (keyed by OG___NNNN IDs; aggregate per OG)
    sig_set, sig_details = collect_signalp()
    tm_counts = collect_tm()
    og_sigcount = Counter()
    og_tmany = Counter()
    og_tmtotal = defaultdict(list)
    for sid in sig_set:
        og = og_of(sid)
        if og:
            og_sigcount[og] += 1
    for sid, n in tm_counts.items():
        og = og_of(sid)
        if og is None:
            continue
        og_tmtotal[og].append(n)
        if n > 0:
            og_tmany[og] += 1

    # Build summary
    summary = []
    for og in ogs:
        total_seq = len(og_seqs[og])
        species = og_species[og]

        # Best BLAST hit by bitscore
        blast_rows = og_blast[og]
        best_blast = max(blast_rows, key=lambda r: r["bitscore"]) if blast_rows else None
        # Most frequent Swiss-Prot subject (by # queries hitting it)
        subj_counter = Counter()
        for r in blast_rows:
            subj_counter[(r["s"], r["desc"])] += 1
        top_subjects = subj_counter.most_common(3)
        blast_cov = len(og_blast_queries[og]) / total_seq if total_seq else 0

        # PFAM top domains
        pfam_rows = og_pfam[og]
        pfam_counter = Counter()
        pfam_desc = {}
        pfam_evalues = defaultdict(list)
        for r in pfam_rows:
            key = (r["target"], r["acc"])
            pfam_counter[key] += 1
            pfam_desc[key] = r["desc"]
            pfam_evalues[key].append(r["evalue"])
        top_pfam = pfam_counter.most_common(5)
        pfam_cov = len(og_pfam_queries[og]) / total_seq if total_seq else 0

        tm_any = og_tmany[og]
        tm_seqs_scored = len(og_tmtotal[og])
        sig_n = og_sigcount[og]

        summary.append({
            "og": og,
            "n": total_seq,
            "species": dict(species),
            "blast_cov": blast_cov,
            "best_blast": best_blast,
            "top_subjects": top_subjects,
            "n_blast_hits": len(blast_rows),
            "pfam_cov": pfam_cov,
            "top_pfam": [(k[0], k[1], pfam_counter[k], pfam_desc[k],
                          min(pfam_evalues[k])) for k in [k for k, _ in top_pfam]],
            "n_pfam_hits": len(pfam_rows),
            "sig_n": sig_n,
            "sig_frac": sig_n / total_seq if total_seq else 0,
            "tm_any": tm_any,
            "tm_scored": tm_seqs_scored,
            "tm_frac": tm_any / tm_seqs_scored if tm_seqs_scored else 0,
        })

    # Write TSV + JSON
    with open(ROOT / "annotation_summary.tsv", "w") as fh:
        fh.write("OG\tn_seqs\tn_species\tspecies\t"
                 "blast_cov\tbest_blast_subject\tbest_blast_desc\tbest_blast_pident\tbest_blast_bits\tbest_blast_evalue\t"
                 "top_blast_subjects\t"
                 "pfam_cov\ttop_pfam\t"
                 "sig_n\tsig_frac\ttm_any\ttm_scored\ttm_frac\n")
        for s in summary:
            best = s["best_blast"]
            best_subj = best["s"] if best else ""
            best_desc = best["desc"].split(" OS=")[0] if best else ""
            best_pid = f"{best['pident']:.1f}" if best else ""
            best_bits = f"{best['bitscore']:.1f}" if best else ""
            best_ev = f"{best['evalue']:.1e}" if best else ""
            top_subj = "; ".join(
                f"{sub.split('|')[-1]} x{n}" for (sub, desc), n in s["top_subjects"]
            )
            top_pfam = "; ".join(
                f"{name}[{acc}] x{n} ({desc}, E={ev:.1e})"
                for name, acc, n, desc, ev in s["top_pfam"]
            )
            species_str = ",".join(f"{k}:{v}" for k, v in sorted(s["species"].items()))
            fh.write(
                f"{s['og']}\t{s['n']}\t{len(s['species'])}\t{species_str}\t"
                f"{s['blast_cov']:.2f}\t{best_subj}\t{best_desc}\t{best_pid}\t{best_bits}\t{best_ev}\t"
                f"{top_subj}\t"
                f"{s['pfam_cov']:.2f}\t{top_pfam}\t"
                f"{s['sig_n']}\t{s['sig_frac']:.2f}\t{s['tm_any']}\t{s['tm_scored']}\t{s['tm_frac']:.2f}\n"
            )

    with open(ROOT / "annotation_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1, default=str)

    print(f"Wrote annotation_summary.tsv and .json for {len(summary)} orthogroups")


if __name__ == "__main__":
    main()
