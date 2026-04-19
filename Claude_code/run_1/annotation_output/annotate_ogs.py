#!/usr/bin/env python3
"""Annotate orthogroups (OGs) by combining BLAST, PFAM, DeepTMHMM and SignalP results.

Outputs:
  annotation_output/og_annotations.tsv  - per-OG summary table
  annotation_output/og_annotations.md   - human-readable report
  annotation_output/og_sequence_table.tsv - per-sequence details
"""
from __future__ import annotations
import json
import os
import re
import glob
from collections import Counter, defaultdict

ROOT = "/Users/xiaoyu/workspace/Orthogroups.calcifying_loose_fastas"
OUT = os.path.join(ROOT, "annotation_output")
os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------------------
# 1. Parse FASTA files - build ordered list of sequence IDs per OG
# ---------------------------------------------------------------------------
def parse_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


og_files = sorted(glob.glob(os.path.join(ROOT, "OG*.fa")))
og_ids = [os.path.basename(p).replace(".fa", "") for p in og_files]
og_to_seqs = {og: parse_fasta_ids(p) for og, p in zip(og_ids, og_files)}

# Build mapping: short_id "OGxxx__NNNNNN" -> original fasta id
short_to_orig = {}
orig_to_short = {}
for og, seqs in og_to_seqs.items():
    for i, sid in enumerate(seqs, 1):
        short = f"{og}__{i:06d}"
        short_to_orig[short] = sid
        orig_to_short[sid] = short


# ---------------------------------------------------------------------------
# 2. Parse BLAST output
# ---------------------------------------------------------------------------
# columns: query, subject_id, subject_full_title, pident, length, evalue, bitscore
def parse_blast(path):
    """Return list of dicts for BLAST hits."""
    hits = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            try:
                pid = float(parts[3])
                evalue = float(parts[5])
                bit = float(parts[6])
            except ValueError:
                continue
            hits.append(
                {
                    "query": parts[0],
                    "subject_id": parts[1],
                    "subject_title": parts[2],
                    "pident": pid,
                    "length": parts[4],
                    "evalue": evalue,
                    "bitscore": bit,
                }
            )
    return hits


og_to_blast = {}
for og in og_ids:
    p = os.path.join(ROOT, "blast_out", f"{og}.blast.tsv")
    og_to_blast[og] = parse_blast(p) if os.path.exists(p) else []


# ---------------------------------------------------------------------------
# 3. Parse PFAM tbl output
# ---------------------------------------------------------------------------
def parse_pfam(path):
    """Return list of dicts for Pfam domain hits."""
    hits = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            # tokens separated by whitespace; description = rest after col 18
            parts = line.rstrip("\n").split(None, 18)
            if len(parts) < 19:
                continue
            try:
                evalue = float(parts[4])
                score = float(parts[5])
            except ValueError:
                continue
            hits.append(
                {
                    "pfam_name": parts[0],
                    "pfam_acc": parts[1],
                    "query": parts[2],
                    "evalue": evalue,
                    "score": score,
                    "description": parts[18].strip(),
                }
            )
    return hits


og_to_pfam = {}
for og in og_ids:
    p = os.path.join(ROOT, "hmm_out", f"{og}.pfam.tbl")
    og_to_pfam[og] = parse_pfam(p) if os.path.exists(p) else []


# ---------------------------------------------------------------------------
# 4. Parse DeepTMHMM 3-line / gff outputs
# ---------------------------------------------------------------------------
TM_STATES = {"M", "B"}  # M = alpha TM helix, B = beta barrel strand


def parse_deeptmhmm_3line(paths):
    """Return dict: short_id -> (label, num_TMRs, length, has_SP)."""
    out = {}
    for path in paths:
        with open(path) as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                if not h.startswith(">"):
                    continue
                header = h[1:].strip()
                # header = "OG0000049__000001 | GLOB"
                sid, _, label = header.partition("|")
                sid = sid.strip()
                label = label.strip()
                _seq = fh.readline().rstrip("\n")
                topo = fh.readline().rstrip("\n")
                tm_count = 0
                # count non-overlapping runs of TM states
                in_tm = False
                for c in topo:
                    if c in TM_STATES:
                        if not in_tm:
                            tm_count += 1
                            in_tm = True
                    else:
                        in_tm = False
                has_sp = "S" in topo
                out[sid] = {
                    "label": label,
                    "tmr_count": tm_count,
                    "length": len(topo),
                    "has_sp": has_sp,
                }
    return out


tmhmm_paths = sorted(
    glob.glob(os.path.join(ROOT, "deeptmhmm_results", "deeptmhmm_batch_*.predicted_topologies.3line"))
)
short_to_tmhmm = parse_deeptmhmm_3line(tmhmm_paths)


# ---------------------------------------------------------------------------
# 5. Parse SignalP JSONs
# ---------------------------------------------------------------------------
short_to_signalp = {}
for jp in sorted(glob.glob(os.path.join(ROOT, "signalp_results", "signalp_batch_*.output.json"))):
    with open(jp) as fh:
        d = json.load(fh)
    for sid, v in d["SEQUENCES"].items():
        pred = v.get("Prediction", "Other")
        cs = v.get("CS_pos", "")
        like = v.get("Likelihood", [])
        short_to_signalp[sid] = {
            "prediction": pred,
            "cs_pos": cs,
            "likelihood_sp": like[1] if len(like) > 1 else 0.0,
        }


# ---------------------------------------------------------------------------
# 6. Calcification relevance knowledge base
# ---------------------------------------------------------------------------
# Keywords / terms commonly associated with calcification processes
# (biomineralization, ion transport, carbonic anhydrase, EF-hand / Ca2+ binding,
#  acidic / matrix proteins, galaxins, pearlins, collagen, cadherins, etc.)
CALC_KEYWORDS = [
    # CaCO3 chemistry / carbon
    ("carbonic anhydrase", "CACA", "carbonate chemistry"),
    ("bicarbonate", "HCO3", "carbonate chemistry"),
    # Ca2+ transport / binding
    ("calcium", "Ca2+", "Ca2+ handling"),
    ("ca2+", "Ca2+", "Ca2+ handling"),
    ("EF-hand", "EF-hand", "Ca2+ binding"),
    ("EF hand", "EF-hand", "Ca2+ binding"),
    ("calmodulin", "CaM", "Ca2+ signaling"),
    ("annexin", "ANX", "Ca2+-phospholipid binding"),
    ("calreticulin", "CRT", "Ca2+ storage"),
    ("calponin", "CNN", "Ca2+ binding"),
    ("sodium/calcium", "NCX", "Ca2+ transport"),
    # transporters relevant to ion flux
    ("SLC4", "SLC4", "bicarbonate transporter"),
    ("SLC26", "SLC26", "anion exchanger"),
    ("V-ATPase", "V-ATPase", "H+ pumping"),
    ("P-type ATPase", "P-ATPase", "ion pump"),
    ("Ca-ATPase", "SERCA/PMCA", "Ca2+ pump"),
    # Biomineralization / matrix
    ("galaxin", "galaxin", "coral matrix protein"),
    ("pearlin", "pearlin", "mollusc nacre protein"),
    ("lustrin", "lustrin", "nacre matrix"),
    ("perlucin", "perlucin", "C-type lectin nacre"),
    ("perlustrin", "perlustrin", "nacre matrix"),
    ("nacrein", "nacrein", "nacre CA-like"),
    ("N16/pearlin", "pearlin", "nacre matrix"),
    ("mucoperlin", "mucoperlin", "nacre matrix"),
    ("calcarin", "calcarin", "foraminifera matrix"),
    ("sm50", "SM50", "sea urchin spicule matrix"),
    ("sm30", "SM30", "sea urchin spicule matrix"),
    ("msp130", "MSP130", "sea urchin mesenchyme"),
    ("c-type lectin", "CTL", "matrix glycoprotein"),
    ("mammalian ependymin", "ependymin", "mineralization related"),
    # Biomineralization structural extracellular proteins
    ("collagen", "COL", "ECM/scaffold"),
    ("fibrillin", "fibrillin", "ECM/microfibril"),
    ("cadherin", "CDH", "cell adhesion"),
    ("chitin", "chitin", "chitinous organic scaffold"),
    ("silicatein", "silicatein", "silica biomineralization"),
    ("laminin", "laminin", "basement membrane ECM"),
    ("fibronectin", "FN", "ECM"),
    ("vitelline", "vitelline", "egg envelope matrix"),
    # Acidic / repetitive mineralization motifs
    ("acidic", "acidic", "acidic matrix protein"),
    ("aspartate rich", "Asp-rich", "acidic matrix"),
    # pH / proton
    ("proton", "H+", "pH regulation"),
    ("chloride", "Cl-", "anion transport"),
    # secretion / vesicle
    ("von Willebrand", "vWF", "extracellular protein"),
    ("sushi", "sushi", "complement/adhesion"),
    # coccolithophore / matrix
    ("adhesive plaque", "adhesive-plaque", "byssal/adhesive matrix protein"),
    ("arylsulfatase", "arylsulfatase", "GAG/sulfated polysaccharide processing"),
    ("sulfotransferase", "sulfotransferase", "sulfation of ECM/GAGs"),
    ("heparan", "heparan", "HS proteoglycan"),
    ("proteoglycan", "proteoglycan", "sulfated ECM"),
    ("coccolith", "coccolith", "coccolithophore shell matrix"),
    ("piccolo", "piccolo", "Ca2+-sensor/active-zone protein"),
    ("shell matrix", "shell-matrix", "shell biomineralization"),
    ("silicalemma", "silicalemma", "silica biomineralization"),
    ("MAP", "mineralization-assoc", "mineralization-associated protein"),
    ("mucin", "mucin", "secreted glycoprotein"),
    ("uromodulin", "uromodulin", "Tamm-Horsfall/ECM"),
]

# PFAM accessions/names commonly associated with biomineralization / Ca handling
CALC_PFAMS = {
    # carbonic anhydrase
    "Carb_anhydrase": "carbonic anhydrase (CA activity)",
    "alpha_CA": "alpha-carbonic anhydrase",
    "Pro_CA": "prokaryotic carbonic anhydrase",
    # Ca2+ binding
    "EF-hand_1": "EF-hand Ca2+-binding domain",
    "EF-hand_2": "EF-hand Ca2+-binding domain",
    "EF-hand_4": "EF-hand Ca2+-binding domain",
    "EF-hand_5": "EF-hand Ca2+-binding domain",
    "EF-hand_6": "EF-hand Ca2+-binding domain",
    "EF-hand_7": "EF-hand Ca2+-binding domain",
    "EF-hand_8": "EF-hand Ca2+-binding domain",
    "EF-hand_9": "EF-hand Ca2+-binding domain",
    "EFhand_Ca_insen": "EF-hand (Ca2+-insensitive)",
    "SPARC_Ca_bdg": "SPARC extracellular Ca2+-binding",
    "Annexin": "annexin (Ca2+-phospholipid)",
    "Calreticulin": "calreticulin (ER Ca2+ storage)",
    "Calpain_III": "calpain (Ca2+ protease)",
    "CaMKII_AD": "CaMKII assoc. domain",
    # Ion transporters
    "Na_Ca_ex": "Na+/Ca2+ exchanger",
    "HCO3_cotransp": "bicarbonate co-transporter",
    "Band_3_cyto": "band3/SLC4",
    "Sulfate_transp": "SLC26 anion exchanger",
    "ATP_Ca_trans_N": "PMCA/SERCA Ca2+ ATPase",
    "Cation_ATPase": "cation ATPase",
    "Cation_ATPase_N": "cation ATPase (P-type)",
    "V_ATPase_I": "V-type ATPase",
    # C-type lectin / nacre-matrix-like
    "Lectin_C": "C-type lectin (matrix glycoprotein)",
    "ricin_B_lectin": "ricin B lectin",
    # Biomineralization-associated proteins
    "Collagen": "collagen triple helix",
    "VWA": "von Willebrand A (ECM/matrix)",
    "VWC": "von Willebrand C (ECM)",
    "VWD": "von Willebrand D",
    "Chitin_bind_1": "chitin-binding peritrophin-A",
    "Chitin_bind_2": "chitin-binding type-2",
    "Chitin_bind_3": "chitin-binding type-3",
    "CBM_14": "chitin-binding CBM14",
    "Cadherin": "cadherin repeat (Ca2+-dep adhesion)",
    "Cadherin_2": "cadherin repeat",
    "Cadherin_3": "cadherin repeat",
    "Laminin_G_1": "laminin G (ECM)",
    "Laminin_G_2": "laminin G",
    "Laminin_G_3": "laminin G",
    "Fibronectin_3": "fibronectin type III (ECM)",
    "Fibronectin_2": "fibronectin type II",
    "Fibronectin_1": "fibronectin type I",
    "fn3": "fibronectin type III",
    "Sushi": "sushi/SCR (complement/adhesion)",
    "SRCR": "scavenger receptor Cys-rich",
    "EGF": "EGF-like (ECM, often Ca2+)",
    "EGF_CA": "Ca2+-binding EGF-like",
    "Thrombospondin": "thrombospondin (ECM)",
    "TSP_1": "thrombospondin type-1 repeat",
    # galaxin-like (coral)
    "Galaxin": "galaxin (coral skeletal matrix)",
    "Pif": "Pif/nacre matrix",
}


def classify_calcification(
    func_text: str,
    pfam_names: list[str],
    pct_sp: float,
    pct_tm: float,
) -> tuple[str, str]:
    """Return (tier, rationale)."""
    text = func_text.lower()
    reasons = []
    # Pfam-driven evidence
    for pname in pfam_names:
        if pname in CALC_PFAMS:
            reasons.append(f"PFAM:{pname} ({CALC_PFAMS[pname]})")
    # keyword evidence (BLAST titles + function text)
    for kw, tag, note in CALC_KEYWORDS:
        if kw.lower() in text:
            reasons.append(f"keyword:{tag} ({note})")
    # secretory / membrane topology hints (supporting, not primary)
    topology_hints = []
    if pct_sp >= 50:
        topology_hints.append(f"predominantly secreted ({pct_sp:.0f}% SP)")
    elif pct_sp >= 20:
        topology_hints.append(f"partially secreted ({pct_sp:.0f}% SP)")
    if pct_tm >= 50:
        topology_hints.append(f"membrane-bound ({pct_tm:.0f}% TM)")
    elif pct_tm >= 20:
        topology_hints.append(f"partial TM ({pct_tm:.0f}% TM)")

    # dedupe preserving order
    seen = set()
    deduped = []
    for r in reasons:
        if r not in seen:
            seen.add(r)
            deduped.append(r)

    strong_markers = [
        "carbonic anhydrase",
        "ef-hand",
        "annexin",
        "calreticulin",
        "pmca",
        "serca",
        "slc4",
        "slc26",
        "na+/ca2+",
        "na_ca_ex",
        "galaxin",
        "nacrein",
        "pearlin",
        "lustrin",
        "perlucin",
        "perlustrin",
        "sm50",
        "sm30",
        "msp130",
        "pif",
        "vwa",
        "vwc",
        "vwd",
        "v-atpase",
        "p-atpase",
        "chitin",
        "collagen",
        "laminin",
        "cadherin",
        "egf_ca",
        "sparc",
        "bicarbonate",
    ]
    medium_markers = [
        "egf",
        "ca2+",
        "calcium",
        "fibronectin",
        "sushi",
        "thrombospondin",
        "apolipoprotein",
        "sulfatase",
        "sulfotransferase",
        "arylsulfatase",
        "heparan",
        "proteoglycan",
        "trypsin",
        "acidic",
        "piccolo",
        "adhesive",
        "mucin",
    ]
    joined = " ".join(deduped).lower()
    tier = "low"
    if any(m in joined for m in strong_markers):
        tier = "high"
    elif any(m in joined for m in medium_markers):
        tier = "medium"
    elif not deduped:
        # no direct evidence but secretion + no TM = candidate secreted matrix protein
        if pct_sp >= 50 and pct_tm < 20:
            tier = "medium"
            deduped.append("secreted protein (candidate extracellular matrix)")

    rationale = "; ".join(deduped) if deduped else "no direct evidence"
    if topology_hints:
        rationale += " | topology: " + ", ".join(topology_hints)
    return tier, rationale


# ---------------------------------------------------------------------------
# 7. Aggregate per OG & build function annotation
# ---------------------------------------------------------------------------
def summarize_titles(titles: list[str]) -> str:
    """Clean and rank blast subject titles."""
    cleaned = []
    for t in titles:
        t = re.sub(r"^sp\|[^|]+\|[^ ]+ ", "", t)  # strip leading sp|.. accession tag
        t = re.sub(r" OS=.*$", "", t)  # drop taxonomy etc.
        t = t.strip()
        cleaned.append(t)
    c = Counter(cleaned)
    return "; ".join(f"{name} (x{cnt})" for name, cnt in c.most_common(3))


def build_function_text(og: str, blast_hits: list, pfam_hits: list) -> tuple[str, str, list[str]]:
    """Return (short_function, long_function, pfam_name_list)."""
    pfam_names = [h["pfam_name"] for h in pfam_hits]
    pfam_counter = Counter(pfam_names)
    pfam_summary = "; ".join(f"{n} (x{c})" for n, c in pfam_counter.most_common(5))

    blast_titles = [h["subject_title"] for h in blast_hits]
    blast_summary = summarize_titles(blast_titles)

    if pfam_summary and blast_summary:
        short = f"{pfam_summary.split(';')[0].strip()} / {blast_summary.split(';')[0].strip()}"
    elif pfam_summary:
        short = pfam_summary.split(";")[0].strip()
    elif blast_summary:
        short = blast_summary.split(";")[0].strip()
    else:
        short = "no detectable homology or domain"

    long_parts = []
    if pfam_summary:
        long_parts.append(f"PFAM: {pfam_summary}")
    if blast_summary:
        long_parts.append(f"BLAST: {blast_summary}")
    long_text = " | ".join(long_parts) if long_parts else "no detectable homology or domain"
    return short, long_text, list(set(pfam_names))


rows = []  # per-OG summary rows
seq_rows = []  # per-sequence detail rows
for og in og_ids:
    seqs = og_to_seqs[og]
    n_seq = len(seqs)
    # per-sequence lookups
    n_sp = 0
    n_tm = 0
    tm_counts = []
    lengths = []
    for i, orig in enumerate(seqs, 1):
        short = f"{og}__{i:06d}"
        sp = short_to_signalp.get(short, {})
        tm = short_to_tmhmm.get(short, {})
        if sp.get("prediction") and sp["prediction"] != "Other":
            n_sp += 1
        c = tm.get("tmr_count", 0)
        tm_counts.append(c)
        if c > 0:
            n_tm += 1
        if tm.get("length"):
            lengths.append(tm["length"])

        seq_rows.append(
            {
                "og": og,
                "short_id": short,
                "orig_id": orig,
                "length": tm.get("length", ""),
                "tmr_count": tm.get("tmr_count", ""),
                "tmhmm_label": tm.get("label", ""),
                "signalp_prediction": sp.get("prediction", ""),
                "signalp_cs_pos": sp.get("cs_pos", ""),
                "signalp_sp_likelihood": f"{sp.get('likelihood_sp', 0):.3f}" if sp else "",
            }
        )

    pct_sp = 100.0 * n_sp / n_seq if n_seq else 0
    pct_tm = 100.0 * n_tm / n_seq if n_seq else 0
    mean_len = sum(lengths) / len(lengths) if lengths else 0
    mean_tm = sum(tm_counts) / len(tm_counts) if tm_counts else 0

    blast = og_to_blast[og]
    pfam = og_to_pfam[og]
    func_short, func_long, pfam_name_list = build_function_text(og, blast, pfam)

    tier, rationale = classify_calcification(func_long, pfam_name_list, pct_sp, pct_tm)

    rows.append(
        {
            "og": og,
            "n_sequences": n_seq,
            "mean_length_aa": f"{mean_len:.0f}",
            "pct_signal_peptide": f"{pct_sp:.1f}",
            "pct_with_TM": f"{pct_tm:.1f}",
            "mean_TMRs": f"{mean_tm:.2f}",
            "n_pfam_hits": len(pfam),
            "n_blast_hits": len(blast),
            "top_pfam": "; ".join(
                f"{n}({c})" for n, c in Counter([h['pfam_name'] for h in pfam]).most_common(5)
            ),
            "top_blast": summarize_titles([h["subject_title"] for h in blast]),
            "function_short": func_short,
            "function_long": func_long,
            "calcification_relevance": tier,
            "calcification_rationale": rationale,
        }
    )


# ---------------------------------------------------------------------------
# 8. Write TSV outputs
# ---------------------------------------------------------------------------
summary_cols = [
    "og",
    "n_sequences",
    "mean_length_aa",
    "pct_signal_peptide",
    "pct_with_TM",
    "mean_TMRs",
    "n_pfam_hits",
    "n_blast_hits",
    "top_pfam",
    "top_blast",
    "function_short",
    "function_long",
    "calcification_relevance",
    "calcification_rationale",
]
with open(os.path.join(OUT, "og_annotations.tsv"), "w") as fh:
    fh.write("\t".join(summary_cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in summary_cols) + "\n")

seq_cols = [
    "og",
    "short_id",
    "orig_id",
    "length",
    "tmr_count",
    "tmhmm_label",
    "signalp_prediction",
    "signalp_cs_pos",
    "signalp_sp_likelihood",
]
with open(os.path.join(OUT, "og_sequence_table.tsv"), "w") as fh:
    fh.write("\t".join(seq_cols) + "\n")
    for r in seq_rows:
        fh.write("\t".join(str(r[c]) for c in seq_cols) + "\n")


# ---------------------------------------------------------------------------
# 9. Write markdown report
# ---------------------------------------------------------------------------
tier_order = {"high": 0, "medium": 1, "low": 2}
rows_sorted = sorted(rows, key=lambda r: (tier_order[r["calcification_relevance"]], r["og"]))

# count species contributions across all OGs
species_counter = Counter()
for og, seqs in og_to_seqs.items():
    for s in seqs:
        parts = s.split("|")
        if len(parts) >= 2:
            species_counter[parts[1]] += 1

with open(os.path.join(OUT, "og_annotations.md"), "w") as fh:
    fh.write("# Orthogroup annotation & calcification relevance\n\n")
    fh.write(
        "Annotations generated from BLASTp (SwissProt), HMMER/Pfam (Pfam-A), "
        "DeepTMHMM topology, and SignalP 6.0 secretion predictions for the "
        "`Orthogroups.calcifying_loose_fastas` set.\n\n"
    )
    fh.write(
        "## Dataset\n\n"
        f"- {len(rows)} orthogroups\n"
        f"- {sum(len(v) for v in og_to_seqs.values())} total sequences across "
        f"{len(species_counter)} source proteomes\n\n"
        "Source proteomes (sequence count) - all haptophytes/coccolithophores "
        "and related chromists used as a *calcifying* query set:\n\n"
    )
    for sp_code, n in species_counter.most_common():
        fh.write(f"- `{sp_code}`: {n}\n")
    fh.write("\n")
    tier_counts = Counter(r["calcification_relevance"] for r in rows)
    fh.write(
        f"- **High** calcification relevance: {tier_counts.get('high', 0)}\n"
        f"- **Medium** relevance: {tier_counts.get('medium', 0)}\n"
        f"- **Low / no evidence**: {tier_counts.get('low', 0)}\n\n"
    )
    fh.write("## Relevance tiers\n\n")
    fh.write(
        "* **high** - PFAM or BLAST match to a canonical calcification / Ca²⁺ / "
        "biomineralization marker (carbonic anhydrase, EF-hand, annexin, SLC4/26, "
        "Ca²⁺-ATPase, galaxin, VWA/VWC, collagen, cadherin, chitin-binding, etc.).\n"
        "* **medium** - weaker or indirect evidence (peripheral ECM/adhesion term or "
        "single keyword match).\n"
        "* **low** - no keyword/domain evidence for calcification (still may be "
        "lineage-specific secreted proteins - see SP/TM columns).\n\n"
    )
    fh.write("## Summary table\n\n")
    fh.write("| OG | Tier | Function | SP% | TM% | Rationale |\n")
    fh.write("|---|---|---|---|---|---|\n")
    for r in rows_sorted:
        func = r["function_short"].replace("|", "/")
        rat = r["calcification_rationale"].replace("|", "/")
        fh.write(
            f"| {r['og']} | {r['calcification_relevance']} | {func} | "
            f"{r['pct_signal_peptide']} | {r['pct_with_TM']} | {rat} |\n"
        )
    fh.write("\n## Per-OG annotations (sorted by relevance)\n\n")
    for r in rows_sorted:
        fh.write(f"### {r['og']}  - _{r['calcification_relevance'].upper()}_\n\n")
        fh.write(
            f"- sequences: {r['n_sequences']}; mean length: {r['mean_length_aa']} aa\n"
            f"- signal peptide: {r['pct_signal_peptide']}%; "
            f"\u22651 TM region: {r['pct_with_TM']}% (mean TMRs {r['mean_TMRs']})\n"
            f"- top PFAM domains: {r['top_pfam'] or '-'}\n"
            f"- top BLAST subjects: {r['top_blast'] or '-'}\n"
            f"- **function**: {r['function_short']}\n"
            f"- calcification rationale: {r['calcification_rationale']}\n\n"
        )

print(f"Wrote {len(rows)} OG annotations.")
print("tier counts:", dict(tier_counts))
