#!/usr/bin/env python3
"""Aggregate BLAST, Pfam, SignalP, DeepTMHMM results per orthogroup."""
import os, re, glob, json, csv, collections, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FA_DIR = ROOT
BLAST_DIR = os.path.join(ROOT, "blast_out")
HMM_DIR = os.path.join(ROOT, "hmm_out")
SIGP_DIR = os.path.join(ROOT, "signalp_results")
TMHMM_DIR = os.path.join(ROOT, "deeptmhmm_results")
OUT_DIR = os.path.join(ROOT, "annotation")

def get_ogs():
    ogs = []
    for f in sorted(os.listdir(FA_DIR)):
        if f.endswith(".fa"):
            ogs.append(f[:-3])
    return ogs

def count_seqs(og):
    n = 0
    with open(os.path.join(FA_DIR, og + ".fa")) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n

def parse_blast(og):
    """Return summary of blast hits for an OG."""
    path = os.path.join(BLAST_DIR, og + ".blast.tsv")
    if not os.path.exists(path):
        return None
    hits = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            q, sub_acc, desc, pident, length, evalue, bitscore = p[:7]
            try:
                e = float(evalue); b = float(bitscore); pid = float(pident); L = int(length)
            except ValueError:
                continue
            # Strip "sp|XX|YY " from description
            clean = re.sub(r"^sp\|[^|]+\|[^ ]+\s+", "", desc)
            # Gene name + organism
            m = re.match(r"(.+?)\s+OS=(.+?)\s+OX=", clean)
            prot_desc = m.group(1) if m else clean
            organism = m.group(2) if m else ""
            gn = re.search(r"GN=([^ ]+)", desc)
            gene = gn.group(1) if gn else ""
            hits.append({
                "query": q, "subject": sub_acc, "desc": prot_desc,
                "gene": gene, "organism": organism, "pident": pid,
                "length": L, "evalue": e, "bitscore": b,
            })
    # Summarize: queries with any hit, most-frequent subject descriptions by bitscore-best-per-query
    per_query_best = {}
    for h in hits:
        q = h["query"]
        if q not in per_query_best or h["bitscore"] > per_query_best[q]["bitscore"]:
            per_query_best[q] = h
    # Count dominant annotations (use protein description)
    desc_count = collections.Counter()
    gene_count = collections.Counter()
    best_evalue = float("inf")
    best_hit = None
    for h in per_query_best.values():
        desc_count[h["desc"]] += 1
        if h["gene"]:
            gene_count[h["gene"]] += 1
        if h["evalue"] < best_evalue:
            best_evalue = h["evalue"]; best_hit = h
    return {
        "n_queries_with_hit": len(per_query_best),
        "n_total_hits": len(hits),
        "top_descs": desc_count.most_common(5),
        "top_genes": gene_count.most_common(5),
        "best_hit": best_hit,
    }

def parse_hmm(og):
    path = os.path.join(HMM_DIR, og + ".pfam.tbl")
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split(None, 18)
            if len(p) < 19:
                continue
            target, acc, qname, qacc, evalue = p[0], p[1], p[2], p[3], p[4]
            desc = p[18].rstrip("\n")
            try:
                e = float(evalue)
            except ValueError:
                continue
            if e > 1e-3:  # E-value cutoff
                continue
            rows.append({"pfam": target, "acc": acc, "query": qname, "evalue": e, "desc": desc})
    per_query_best = {}
    domain_count = collections.Counter()
    domain_desc = {}
    for r in rows:
        q = r["query"]
        if q not in per_query_best or r["evalue"] < per_query_best[q]["evalue"]:
            per_query_best[q] = r
    # Count per query which domains they have (any hit)
    for q, r in per_query_best.items():
        domain_count[r["pfam"]] += 1
        domain_desc[r["pfam"]] = r["desc"]
    # Also count across all rows for multi-domain situations
    any_domain_queries = set(r["query"] for r in rows)
    multi_domain_count = collections.Counter()
    for r in rows:
        multi_domain_count[r["pfam"]] += 1
    multi_domain_desc = {r["pfam"]: r["desc"] for r in rows}
    return {
        "n_queries_with_hit": len(any_domain_queries),
        "top_domains": [(d, c, multi_domain_desc.get(d, "")) for d, c in multi_domain_count.most_common(6)],
    }

def parse_signalp():
    """Return dict: seq_id -> prediction (SP probability)."""
    out = {}
    for path in glob.glob(os.path.join(SIGP_DIR, "*_all_results/prediction_results.txt")):
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                sid = p[0]; pred = p[1]
                try:
                    sp_prob = float(p[3])
                except ValueError:
                    sp_prob = 0.0
                # Batch 002 supersedes 001 if both have, but prefer SP+ call
                out[sid] = (pred, sp_prob)
    return out

def parse_tmhmm():
    """Return dict: seq_id -> n_TMRs."""
    out = {}
    for path in sorted(glob.glob(os.path.join(TMHMM_DIR, "*.TMRs.gff3"))):
        cur_id = None
        cur_n = None
        with open(path) as fh:
            for line in fh:
                if line.startswith("# "):
                    m = re.match(r"#\s+(\S+)\s+Number of predicted TMRs:\s+(\d+)", line)
                    if m:
                        out[m.group(1)] = int(m.group(2))
    return out

def main():
    ogs = get_ogs()
    sigp = parse_signalp()
    tmhmm = parse_tmhmm()

    summary_rows = []
    for og in ogs:
        n_seqs = count_seqs(og)
        # Count SP+ & TM+ per OG
        og_sigp_positive = 0
        og_sigp_total = 0
        for sid, (pred, prob) in sigp.items():
            if sid.startswith(og + "__"):
                og_sigp_total += 1
                if pred != "OTHER":
                    og_sigp_positive += 1
        og_tm_positive = 0
        og_tm_counts = []
        og_tm_total = 0
        for sid, n in tmhmm.items():
            if sid.startswith(og + "__"):
                og_tm_total += 1
                og_tm_counts.append(n)
                if n > 0:
                    og_tm_positive += 1
        blast = parse_blast(og)
        hmm = parse_hmm(og)
        summary_rows.append({
            "og": og,
            "n_seqs": n_seqs,
            "blast": blast,
            "hmm": hmm,
            "sigp_pos": og_sigp_positive,
            "sigp_total": og_sigp_total,
            "tm_pos": og_tm_positive,
            "tm_total": og_tm_total,
            "tm_mean": (sum(og_tm_counts) / len(og_tm_counts)) if og_tm_counts else 0.0,
            "tm_max": max(og_tm_counts) if og_tm_counts else 0,
        })

    # Write TSV
    tsv_path = os.path.join(OUT_DIR, "og_annotation_summary.tsv")
    with open(tsv_path, "w") as out:
        out.write("OG\tN_seqs\tBLAST_hits_queries\tTop_BLAST_desc\tBest_BLAST_gene\tBest_BLAST_evalue\tBest_BLAST_organism\tPfam_hits_queries\tTop_Pfam_domains\tSigP_pos/total\tTMR_pos/total\tMean_TMRs\tMax_TMRs\n")
        for r in summary_rows:
            b = r["blast"]; h = r["hmm"]
            top_desc = ""
            best_gene = ""; best_e = ""; best_org = ""
            if b and b["best_hit"]:
                top_desc = "; ".join([f"{d}(n={c})" for d, c in b["top_descs"][:3]])
                bh = b["best_hit"]
                best_gene = bh["gene"]; best_e = f"{bh['evalue']:.1e}"; best_org = bh["organism"]
            top_dom = ""
            if h:
                top_dom = "; ".join([f"{d}(n={c}:{desc[:40]})" for d, c, desc in h["top_domains"][:3]])
            out.write("\t".join([
                r["og"], str(r["n_seqs"]),
                str(b["n_queries_with_hit"]) if b else "0",
                top_desc, best_gene, best_e, best_org,
                str(h["n_queries_with_hit"]) if h else "0",
                top_dom,
                f"{r['sigp_pos']}/{r['sigp_total']}",
                f"{r['tm_pos']}/{r['tm_total']}",
                f"{r['tm_mean']:.2f}",
                str(r["tm_max"]),
            ]) + "\n")
    print("Wrote", tsv_path)

    # Also write a detailed per-OG markdown block for manual annotation
    md_path = os.path.join(OUT_DIR, "og_details.md")
    with open(md_path, "w") as out:
        for r in summary_rows:
            og = r["og"]; b = r["blast"]; h = r["hmm"]
            out.write(f"## {og} (n={r['n_seqs']})\n\n")
            out.write(f"- SignalP positive: {r['sigp_pos']}/{r['sigp_total']}\n")
            out.write(f"- DeepTMHMM with TMRs: {r['tm_pos']}/{r['tm_total']} (mean {r['tm_mean']:.2f}, max {r['tm_max']})\n")
            if b:
                out.write(f"- BLAST: {b['n_queries_with_hit']}/{r['n_seqs']} queries with SwissProt hit\n")
                out.write("  - Top descriptions:\n")
                for d, c in b["top_descs"][:5]:
                    out.write(f"    - {d} (n={c})\n")
                if b["top_genes"]:
                    out.write("  - Top genes: " + ", ".join([f"{g}(n={c})" for g, c in b["top_genes"][:5]]) + "\n")
                if b["best_hit"]:
                    bh = b["best_hit"]
                    out.write(f"  - Best hit: {bh['desc']} | gene={bh['gene']} | {bh['organism']} | E={bh['evalue']:.1e} | pid={bh['pident']:.1f}%\n")
            else:
                out.write(f"- BLAST: no hits file\n")
            if h:
                out.write(f"- Pfam: {h['n_queries_with_hit']}/{r['n_seqs']} queries with domain hit (E<=1e-3)\n")
                for d, c, desc in h["top_domains"][:6]:
                    out.write(f"    - {d} (n={c}) — {desc}\n")
            else:
                out.write(f"- Pfam: no hits file\n")
            out.write("\n")
    print("Wrote", md_path)

if __name__ == "__main__":
    main()
