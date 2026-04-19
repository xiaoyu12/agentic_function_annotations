#!/usr/bin/env python3
"""Annotate orthogroups by aggregating BLAST, Pfam, SignalP, and DeepTMHMM results.

Produces:
- annotations/per_og/<OG>.txt   : detailed per-OG annotation
- annotations/og_summary.tsv    : machine-readable master table
- annotations/og_summary.md     : human-readable report grouped by calcification relevance
"""
from __future__ import annotations

import glob
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path("/Users/xiaoyu/workspace/Orthogroups.calcifying_loose_fastas")
OUT = ROOT / "annotations"
(OUT / "per_og").mkdir(parents=True, exist_ok=True)


# ---------- helpers ----------

def list_ogs() -> list[str]:
    return sorted(p.stem for p in ROOT.glob("OG*.fa"))


def read_fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = None
    seq_parts: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts)))
                name = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq_parts)))
    return records


def composition(seq: str) -> dict[str, float]:
    n = len(seq) or 1
    aa = seq.upper()
    d = aa.count("D")
    e = aa.count("E")
    s = aa.count("S")
    t = aa.count("T")
    c = aa.count("C")
    k = aa.count("K")
    r = aa.count("R")
    p = aa.count("P")
    g = aa.count("G")
    return {
        "length": len(seq),
        "acidic_pct": 100.0 * (d + e) / n,
        "D_pct": 100.0 * d / n,
        "E_pct": 100.0 * e / n,
        "ST_pct": 100.0 * (s + t) / n,
        "basic_pct": 100.0 * (k + r) / n,
        "C_pct": 100.0 * c / n,
        "G_pct": 100.0 * g / n,
        "P_pct": 100.0 * p / n,
    }


def build_og_id_map() -> dict[str, dict[str, str]]:
    """For each OG, map rename-id (OG____000N) -> original protein id."""
    mapping: dict[str, dict[str, str]] = {}
    for og in list_ogs():
        ids = read_fasta_ids(ROOT / f"{og}.fa")
        mapping[og] = {f"{og}__{i+1:06d}": pid for i, pid in enumerate(ids)}
    return mapping


# ---------- BLAST ----------

def parse_blast(path: Path) -> dict[str, list[dict]]:
    """protein_id -> list of hit dicts (top by evalue)."""
    hits: dict[str, list[dict]] = defaultdict(list)
    if not path.exists():
        return hits
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            q, sseqid, stitle, pident, length, evalue, bitscore = parts[:7]
            try:
                e = float(evalue)
                b = float(bitscore)
                p = float(pident)
            except ValueError:
                continue
            hits[q].append({
                "sseqid": sseqid,
                "title": stitle,
                "pident": p,
                "length": int(length) if length.isdigit() else 0,
                "evalue": e,
                "bitscore": b,
            })
    for q in hits:
        hits[q].sort(key=lambda h: (h["evalue"], -h["bitscore"]))
    return hits


# ---------- Pfam ----------

def parse_pfam(path: Path) -> dict[str, list[dict]]:
    """protein_id -> list of domain hits (sorted by evalue)."""
    res: dict[str, list[dict]] = defaultdict(list)
    if not path.exists():
        return res
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split(None, 18)
            if len(parts) < 19:
                continue
            tname, tacc, qname, qacc, e_full, s_full, bias = parts[:7]
            desc = parts[18].rstrip()
            try:
                e = float(e_full)
                s = float(s_full)
            except ValueError:
                continue
            res[qname].append({
                "pfam": tname,
                "accession": tacc,
                "evalue": e,
                "score": s,
                "desc": desc,
            })
    for q in res:
        res[q].sort(key=lambda h: h["evalue"])
    return res


# ---------- SignalP ----------

def parse_signalp() -> dict[str, dict]:
    """rename_id -> {'prediction': str, 'cs_pos': str, 'likelihood': float}."""
    preds: dict[str, dict] = {}
    for p in ROOT.glob("signalp_results/signalp_batch_*.output.json"):
        data = json.loads(p.read_text())
        for name, rec in data.get("SEQUENCES", {}).items():
            preds[name] = {
                "prediction": rec.get("Prediction", "Other"),
                "cs_pos": rec.get("CS_pos", ""),
                "likelihood": max(rec.get("Likelihood", [0]) or [0]),
            }
    return preds


def parse_signalp_gff() -> dict[str, tuple[int, int, float]]:
    """rename_id -> (start, end, score) for predicted signal peptides."""
    out: dict[str, tuple[int, int, float]] = {}
    for p in ROOT.glob("signalp_results/signalp_batch_*.output_all_results/output.gff3"):
        with p.open() as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip().split("\t")
                if len(parts) < 6:
                    continue
                name, _src, feat, start, end, score = parts[:6]
                if feat != "signal_peptide":
                    continue
                try:
                    out[name] = (int(start), int(end), float(score))
                except ValueError:
                    continue
    return out


# ---------- DeepTMHMM ----------

def parse_tmhmm() -> dict[str, dict]:
    """rename_id -> {'n_tmrs': int, 'topology': str, 'length': int}."""
    out: dict[str, dict] = {}
    # topology file has ">name | TYPE" then seq then topology string per protein
    for p in ROOT.glob("deeptmhmm_results/deeptmhmm_batch_*.predicted_topologies.3line"):
        with p.open() as fh:
            block: list[str] = []
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if len(block) == 3:
                        _process_tmhmm_block(block, out)
                    header = line[1:]
                    name = header.split("|")[0].strip()
                    typ = header.split("|")[1].strip() if "|" in header else ""
                    block = [name, typ, ""]
                elif block:
                    # 2nd non-header line is sequence, 3rd is topology
                    if len(block) == 3:
                        # This is the sequence line appended via string
                        block.append(line)
                    else:
                        block.append(line)
            if len(block) >= 5:
                _process_tmhmm_block(block, out)
    # Also use the TMRs gff3 for TMR counts
    for p in ROOT.glob("deeptmhmm_results/deeptmhmm_batch_*.TMRs.gff3"):
        with p.open() as fh:
            current = None
            length = 0
            n_tmrs = 0
            for line in fh:
                if line.startswith("# ") and "Length:" in line:
                    m = re.match(r"#\s+(\S+)\s+Length:\s+(\d+)", line)
                    if m:
                        current = m.group(1)
                        length = int(m.group(2))
                elif line.startswith("# ") and "Number of predicted TMRs" in line:
                    m = re.match(r"#\s+(\S+)\s+Number of predicted TMRs:\s+(\d+)", line)
                    if m and current == m.group(1):
                        n_tmrs = int(m.group(2))
                        rec = out.setdefault(current, {"type": "", "topology": ""})
                        rec["n_tmrs"] = n_tmrs
                        rec["length"] = length
    return out


def _process_tmhmm_block(block: list[str], out: dict[str, dict]):
    # block = [name, type, seq_placeholder_or_seq, seq_or_topo, topo_or_extra]
    # Actually structure appended: [name, type, "", seq, topo]
    if len(block) < 5:
        return
    name, typ, _, _seq, topo = block[:5]
    out[name] = {"type": typ, "topology": topo, "length": len(topo)}


# ---------- calcification scoring ----------

CALC_KEYWORDS = [
    # Keyword patterns target BLAST descriptions or Pfam domain descriptions only.
    # The short accession field is excluded to avoid false positives from UniProt codes like "AE3_CULQU".
    (re.compile(r"carbonic anhydrase", re.I), "carbonic_anhydrase", 4),
    (re.compile(r"calcium[-\s]transporting|calcium[-\s]translocating", re.I), "ca_pump", 4),
    (re.compile(r"\bSERCA\b|\bPMCA\b|plasma membrane calcium", re.I), "ca_pump", 4),
    (re.compile(r"bicarbonate transport|\bSLC4[A-Z]?\d?\b|\bSLC26[A-Z]?\d?\b|anion exchanger \d", re.I), "bicarbonate_anion_exch", 4),
    (re.compile(r"sodium/calcium exchang|\bNCX\d?\b|\bNCKX\d?\b", re.I), "na_ca_exchanger", 4),
    (re.compile(r"calcium channel|voltage[-\s]dependent calcium", re.I), "ca_channel", 3),
    (re.compile(r"\bTRPV\d?\b|\bTRPP\d?\b|polycystin|otopetr", re.I), "ca_channel_trp", 3),
    (re.compile(r"V-?type proton ATPase|V-?ATPase|vacuolar (?:H\+|proton)", re.I), "v_atpase", 2),
    (re.compile(r"calcium[-\s]sensing receptor|\bCaSR\b", re.I), "ca_sensing_receptor", 4),
    (re.compile(r"coccolith|\bGPA\b|\bCruciplacolith\b", re.I), "coccolith_protein", 4),
    (re.compile(r"ef[-\s]hand", re.I), "ef_hand", 2),
    (re.compile(r"calmodulin|calcineurin|centrin\b", re.I), "calmodulin_like", 2),
    (re.compile(r"\bC2 domain|C2-domain", re.I), "c2_domain", 1),
    (re.compile(r"annexin", re.I), "annexin", 3),
    (re.compile(r"galaxin|starmaker|otolin|osteopontin|osteocalcin", re.I), "biomineral_matrix", 4),
    (re.compile(r"aspartic[-\s]rich|glutamic[-\s]rich|acidic.{0,20}(?:region|repeat|rich)", re.I), "acidic_rich", 3),
    (re.compile(r"\bchitin|chitin synth|cellulose synth|glucan synth", re.I), "polysacch_synth", 2),
    (re.compile(r"glycosyltransferase|glycosyl transferase", re.I), "glycosyltransferase", 1),
    (re.compile(r"sulfotransferase|sulphotransferase|arylsulfatase|sulfatase", re.I), "sulfation", 1),
    (re.compile(r"mucin|mucilage", re.I), "mucin_like", 2),
    (re.compile(r"c[-\s]type lectin|c-lectin|mannan-binding", re.I), "lectin", 1),
    (re.compile(r"cadherin", re.I), "cadherin", 1),
    (re.compile(r"alkaline phosphatase", re.I), "alkaline_phosphatase", 1),
    (re.compile(r"magnesium transport|\bCNNM\d?\b|MgtE", re.I), "mg_transport", 2),
    (re.compile(r"phosphate transport|\bNaPi\b|\bPiT\d\b|SLC20", re.I), "pi_transport", 1),
    (re.compile(r"carbonic anhydr|CAH\d|\bALKP\b", re.I), "ca_related", 2),
]


def score_calcification(blast_hits: list[dict], pfam_hits: list[dict]) -> tuple[int, set[str]]:
    tags: set[str] = set()
    score = 0
    for h in blast_hits[:5]:
        # Strip the leading "sp|X|Y" accession so only the human-readable description is searched.
        title = h.get("title") or ""
        title = re.sub(r"^sp\|[^|]+\|\S+\s+", "", title)
        for rx, tag, w in CALC_KEYWORDS:
            if rx.search(title):
                if tag not in tags:
                    score += w
                tags.add(tag)
    for d in pfam_hits:
        text = d.get("desc", "") + " " + d.get("pfam", "")
        for rx, tag, w in CALC_KEYWORDS:
            if rx.search(text):
                if tag not in tags:
                    score += w
                tags.add(tag)
    return score, tags


# ---------- aggregation ----------

def majority(values: list[str]) -> str:
    if not values:
        return ""
    c = Counter(values)
    top, n = c.most_common(1)[0]
    return f"{top} ({n}/{len(values)})"


def top_titles(hits_by_prot: dict[str, list[dict]]) -> list[tuple[str, int]]:
    """Count how often each swissprot title appears as the top hit across members."""
    counter: Counter = Counter()
    for prot, hits in hits_by_prot.items():
        if not hits:
            continue
        title = hits[0].get("title", "").split("OS=")[0].strip()
        # strip leading accession prefix like sp|X|Y
        title = re.sub(r"^sp\|[^|]+\|\S+\s+", "", title)
        counter[title] += 1
    return counter.most_common()


def top_pfams(pf_by_prot: dict[str, list[dict]]) -> list[tuple[str, int, str]]:
    counter: Counter = Counter()
    desc_map: dict[str, str] = {}
    for prot, hits in pf_by_prot.items():
        seen_here: set[str] = set()
        for h in hits:
            key = h["pfam"]
            if key in seen_here:
                continue
            seen_here.add(key)
            counter[key] += 1
            desc_map[key] = h["desc"]
    return [(pf, n, desc_map.get(pf, "")) for pf, n in counter.most_common()]


def main():
    ogs = list_ogs()
    og_map = build_og_id_map()
    sp_preds = parse_signalp()
    sp_gff = parse_signalp_gff()
    tm_preds = parse_tmhmm()

    master_rows: list[dict] = []

    for og in ogs:
        records = read_fasta(ROOT / f"{og}.fa")
        fasta_ids = [r[0] for r in records]
        seqs = {r[0]: r[1] for r in records}
        n_total = len(fasta_ids)

        blast = parse_blast(ROOT / "blast_out" / f"{og}.blast.tsv")
        pfam = parse_pfam(ROOT / "hmm_out" / f"{og}.pfam.tbl")

        # Map rename IDs to originals for SignalP/TMHMM lookup
        rename_to_orig = og_map[og]

        sp_hits = 0
        tm_multi = 0
        tm_single = 0
        tm_any = 0
        cs_positions: list[str] = []
        tm_types: list[str] = []

        for rn, orig in rename_to_orig.items():
            if rn in sp_preds and sp_preds[rn]["prediction"] != "Other":
                sp_hits += 1
                if sp_preds[rn].get("cs_pos"):
                    cs_positions.append(sp_preds[rn]["cs_pos"])
            rec = tm_preds.get(rn)
            if rec:
                n = rec.get("n_tmrs", 0) or 0
                if n > 0:
                    tm_any += 1
                if n == 1:
                    tm_single += 1
                elif n > 1:
                    tm_multi += 1
                if rec.get("type"):
                    tm_types.append(rec["type"])

        prot_with_blast = sum(1 for pid in fasta_ids if pid in blast)
        prot_with_pfam = sum(1 for pid in fasta_ids if pid in pfam)

        titles = top_titles({pid: blast.get(pid, []) for pid in fasta_ids})
        pfams = top_pfams({pid: pfam.get(pid, []) for pid in fasta_ids})

        # Calcification scoring: aggregate across all members
        total_score = 0
        all_tags: set[str] = set()
        for pid in fasta_ids:
            s, t = score_calcification(blast.get(pid, []), pfam.get(pid, []))
            total_score += s
            all_tags.update(t)
        # Normalise: avoid giving huge OGs crazy scores — use 'at least one hit' tags
        calc_score = len(all_tags) * 2 + min(total_score, 50)

        # --- sequence composition / biomineralization candidate flags ---
        comps = [composition(seqs[pid]) for pid in fasta_ids]
        mean_len = sum(c["length"] for c in comps) / max(1, len(comps))
        mean_acidic = sum(c["acidic_pct"] for c in comps) / max(1, len(comps))
        mean_ST = sum(c["ST_pct"] for c in comps) / max(1, len(comps))
        mean_cys = sum(c["C_pct"] for c in comps) / max(1, len(comps))

        # Biomineralization/secretome candidate: signal-peptide + acidic-rich, or uncharacterized + secreted
        bmin_flags: list[str] = []
        pct_sp = 100 * sp_hits / n_total if n_total else 0
        if pct_sp >= 50 and mean_acidic >= 15:
            bmin_flags.append("secreted_acidic")
        elif pct_sp >= 50:
            bmin_flags.append("secretome_candidate")
        if pct_sp >= 30 and prot_with_blast / max(1, n_total) < 0.2:
            bmin_flags.append("novel_secreted")
        if mean_acidic >= 18:
            bmin_flags.append("acidic_rich")
        if mean_cys >= 5:
            bmin_flags.append("cysteine_rich")
        pct_tm_multi = 100 * tm_multi / n_total if n_total else 0
        if pct_tm_multi >= 60 and prot_with_blast / max(1, n_total) < 0.3:
            bmin_flags.append("novel_multipass_membrane")

        if bmin_flags:
            # Boost the calcification score if the OG looks like a biomineralization candidate
            calc_score += 2 * len(bmin_flags)

        # Compose a consensus annotation
        consensus = titles[0][0] if titles else "no annotation"
        consensus_pfam = ", ".join(f"{p} ({n}x: {d})" for p, n, d in pfams[:3])

        # --- write per-OG file ---
        per_path = OUT / "per_og" / f"{og}.txt"
        with per_path.open("w") as fh:
            fh.write(f"Orthogroup: {og}\n")
            fh.write(f"Members: {n_total}\n")
            fh.write(f"With BLAST hit: {prot_with_blast}  |  With Pfam: {prot_with_pfam}\n")
            fh.write(f"SignalP positive: {sp_hits}/{n_total}\n")
            fh.write(f"TMR: any={tm_any}/{n_total}  single-pass={tm_single}  multi-pass={tm_multi}\n")
            if tm_types:
                fh.write(f"TMHMM topology types: {Counter(tm_types).most_common()}\n")
            fh.write("\nTop SwissProt hit(s) across members:\n")
            for title, n in titles[:10]:
                fh.write(f"  [{n}x] {title}\n")
            fh.write("\nPfam domains observed:\n")
            for pf, n, d in pfams[:15]:
                fh.write(f"  [{n}x] {pf} — {d}\n")
            fh.write("\nCalcification tags: " + ", ".join(sorted(all_tags)) + "\n")
            fh.write(f"Biomin flags: {','.join(bmin_flags) or '-'}\n")
            fh.write(f"Mean length: {mean_len:.0f} aa | mean acidic (D+E): {mean_acidic:.1f}% | "
                     f"mean S+T: {mean_ST:.1f}% | mean Cys: {mean_cys:.1f}%\n")
            fh.write(f"Calcification score: {calc_score}\n")

        master_rows.append({
            "OG": og,
            "n_members": n_total,
            "n_blast": prot_with_blast,
            "n_pfam": prot_with_pfam,
            "n_signalp": sp_hits,
            "pct_signalp": round(100 * sp_hits / n_total, 1) if n_total else 0,
            "n_tm_any": tm_any,
            "n_tm_single": tm_single,
            "n_tm_multi": tm_multi,
            "pct_tm": round(100 * tm_any / n_total, 1) if n_total else 0,
            "mean_len": round(mean_len, 0),
            "acidic_pct": round(mean_acidic, 1),
            "cys_pct": round(mean_cys, 2),
            "top_hit": consensus,
            "top_hit_count": titles[0][1] if titles else 0,
            "top_pfams": consensus_pfam,
            "calc_tags": ",".join(sorted(all_tags)) or "-",
            "bmin_flags": ",".join(bmin_flags) or "-",
            "calc_score": calc_score,
        })

    # --- write master TSV ---
    tsv = OUT / "og_summary.tsv"
    with tsv.open("w") as fh:
        cols = ["OG", "n_members", "n_blast", "n_pfam", "n_signalp", "pct_signalp",
                "n_tm_any", "n_tm_single", "n_tm_multi", "pct_tm",
                "mean_len", "acidic_pct", "cys_pct",
                "calc_score", "calc_tags", "bmin_flags",
                "top_hit", "top_hit_count", "top_pfams"]
        fh.write("\t".join(cols) + "\n")
        for row in sorted(master_rows, key=lambda r: (-r["calc_score"], r["OG"])):
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")

    # --- write markdown report ---
    md = OUT / "og_summary.md"
    sorted_rows = sorted(master_rows, key=lambda r: (-r["calc_score"], r["OG"]))
    with md.open("w") as fh:
        fh.write("# Orthogroup annotation summary\n\n")
        fh.write(f"Analysed {len(master_rows)} orthogroups. ")
        fh.write("Calcification relevance scored by keyword matches to BLAST titles and Pfam descriptions ")
        fh.write("(carbonic anhydrases, Ca/HCO3 transporters, EF-hand/annexin, coccolith-associated proteins, etc.).\n\n")

        def bucket(r):
            if r["calc_score"] >= 10:
                return "High"
            if r["calc_score"] >= 5:
                return "Medium"
            if r["calc_score"] >= 2:
                return "Low"
            return "Background"

        for label in ("High", "Medium", "Low", "Background"):
            subset = [r for r in sorted_rows if bucket(r) == label]
            if not subset:
                continue
            fh.write(f"## {label} calcification relevance ({len(subset)})\n\n")
            fh.write("| OG | n | SigP% | TM% | D+E% | Calc tags | Biomin flags | Top hit | Top Pfam |\n")
            fh.write("|---|---|---|---|---|---|---|---|---|\n")
            for r in subset:
                tags = r["calc_tags"]
                top_hit = (r["top_hit"] or "-")[:80].replace("|", "\\|")
                top_pf = (r["top_pfams"] or "-")[:80].replace("|", "\\|")
                fh.write(
                    f"| {r['OG']} | {r['n_members']} | {r['pct_signalp']} | {r['pct_tm']} | "
                    f"{r['acidic_pct']} | {tags} | {r['bmin_flags']} | {top_hit} | {top_pf} |\n"
                )
            fh.write("\n")

    print(f"wrote {tsv}")
    print(f"wrote {md}")
    print(f"wrote {len(master_rows)} per-OG files under {OUT/'per_og'}")


if __name__ == "__main__":
    main()
