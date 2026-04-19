#!/usr/bin/env python3
"""Annotate orthogroups from blast, PFAM, DeepTMHMM, SignalP outputs
and rate their relevance to biomineralization / calcification."""
import os, json, re, csv
from collections import Counter, defaultdict

ROOT = '/Users/xiaoyu/workspace/Orthogroups.calcifying_loose_fastas'

# ---------------------------------------------------------------------------
# 1. Build ID map: OG__NNNN -> original fasta header
# ---------------------------------------------------------------------------
og_seqs = {}
ogs = sorted(f[:-3] for f in os.listdir(ROOT) if f.endswith('.fa'))
for og in ogs:
    with open(os.path.join(ROOT, og + '.fa')) as fh:
        og_seqs[og] = [l[1:].split()[0] for l in fh if l.startswith('>')]

# ---------------------------------------------------------------------------
# 2. Parse BLAST outputs
# ---------------------------------------------------------------------------
og_blast = defaultdict(list)
for og in ogs:
    path = os.path.join(ROOT, 'blast_out', f'{og}.blast.tsv')
    if not os.path.exists(path):
        continue
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue
            qid, sid, sdesc, pid, aln, ev, bits = parts[:7]
            og_blast[og].append({
                'qid': qid, 'sid': sid, 'desc': sdesc,
                'pident': float(pid), 'len': int(aln),
                'evalue': float(ev), 'bits': float(bits)
            })

# ---------------------------------------------------------------------------
# 3. Parse PFAM tblout
# ---------------------------------------------------------------------------
og_pfam = defaultdict(list)
for og in ogs:
    path = os.path.join(ROOT, 'hmm_out', f'{og}.pfam.tbl')
    if not os.path.exists(path):
        continue
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            toks = line.rstrip('\n').split()
            if len(toks) < 19:
                continue
            og_pfam[og].append({
                'pfam': toks[0], 'acc': toks[1], 'qid': toks[2],
                'evalue': float(toks[4]), 'score': float(toks[5]),
                'desc': ' '.join(toks[18:])
            })

# ---------------------------------------------------------------------------
# 4. DeepTMHMM topology labels per sequence
# ---------------------------------------------------------------------------
og_tm = defaultdict(Counter)
for i in range(1, 8):
    path = os.path.join(ROOT, 'deeptmhmm_results',
                        f'deeptmhmm_batch_{i:03d}.predicted_topologies.3line')
    if not os.path.exists(path):
        continue
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                hdr = line[1:].strip()
                parts = [p.strip() for p in hdr.split('|')]
                og_base = parts[0].split('__')[0]
                label = parts[1] if len(parts) > 1 else 'NA'
                og_tm[og_base][label] += 1

# ---------------------------------------------------------------------------
# 5. SignalP predictions
# ---------------------------------------------------------------------------
og_sp = defaultdict(Counter)
for i in range(1, 3):
    path = os.path.join(ROOT, 'signalp_results', f'signalp_batch_{i:03d}.output.json')
    if not os.path.exists(path):
        continue
    with open(path) as fh:
        data = json.load(fh)
    for ognum, rec in data['SEQUENCES'].items():
        og_base = ognum.split('__')[0]
        og_sp[og_base][rec['Prediction']] += 1

# ---------------------------------------------------------------------------
# 6. Helpers
# ---------------------------------------------------------------------------
def short_desc(d):
    d = re.sub(r'\s+OS=.*$', '', d)
    d = re.sub(r'^sp\|[^|]+\|\S+\s+', '', d)
    return d.strip()

def consensus_blast(hits, k=3):
    best = {}
    for h in hits:
        s = h['sid']
        if s not in best or h['bits'] > best[s]['bits']:
            best[s] = h
    return sorted(best.values(), key=lambda h: -h['bits'])[:k]

def consensus_pfam(hits, k=6):
    by = defaultdict(lambda: {'n': 0, 'best_ev': 1.0, 'desc': '', 'acc': ''})
    seen = set()
    for h in hits:
        key = (h['pfam'], h['qid'])
        if key in seen: continue
        seen.add(key)
        r = by[h['pfam']]
        r['n'] += 1
        r['desc'] = h['desc']; r['acc'] = h['acc']
        r['best_ev'] = min(r['best_ev'], h['evalue'])
    return sorted(by.items(), key=lambda kv: (-kv[1]['n'], kv[1]['best_ev']))[:k]

# ---------------------------------------------------------------------------
# 7. Curated per-OG functional annotation and calcification-relevance
# ---------------------------------------------------------------------------
# relevance tiers: 'strong', 'moderate', 'weak', 'unlikely'
# Based on PFAM/BLAST evidence + known biomineralization literature
CURATED = {
    # ---- no clear calcification link ----
    'OG0000049': ('FKBP-15-like peptidyl-prolyl isomerase / coiled-coil scaffolding protein (large 502-seq family, mostly in Callep/Pleelo)',
                  'unlikely',
                  'FKBP/coiled-coil functions (cytoskeleton, protein folding) — no established calcification link; note 495/502 globular.'),
    'OG0001332': ('Hint-domain / Hedgehog-like auto-processing proteins (Warthog/Desert-hedgehog blast); 73/91 with TM',
                  'moderate',
                  'Hedgehog signalling autoprocessing domain. In metazoans hedgehog regulates osteoblast differentiation/endochondral ossification; in protists Hint proteins can process adhesive peptides. Possible role in cell-matrix signalling relevant to biomineralisation; high TM fraction suggests membrane-anchored processors.'),
    'OG0001976': ('Pentapeptide-repeat / decapeptide-repeat proteins (72 seqs with PF00805); blast hits to cyanobacterial slr proteins',
                  'moderate',
                  'Pentapeptide repeats form right-handed quadrilateral β-helices; a subfamily (RfrA-like) binds Ca2+/regulates inorganic carbon uptake in cyanobacteria and has been reported in algal CCM contexts. Relevance via carbon concentration / Ca2+ handling plausible but indirect.'),
    'OG0002887': ('PIF1-like DNA helicase (nuclear/mitochondrial genome maintenance)',
                  'unlikely',
                  'DNA helicase — housekeeping role, no calcification link.'),
    'OG0003717': ('Uncharacterised (no blast/PFAM hit); mostly globular, 4 TM',
                  'weak',
                  'No functional evidence; lineage-restricted hypothetical protein. Cannot be assessed.'),
    'OG0009246': ('Asparaginyl-tRNA synthetase (PF00152) but blast picks up Adhesive-plaque-matrix / Collagen-alpha-5(VI); 11/31 have signal peptide',
                  'moderate',
                  'BLAST to mussel adhesive-plaque-matrix and type-VI collagen suggests a secreted matrix-like protein with an embedded aaRS domain. High signal-peptide fraction is consistent with an extracellular matrix role relevant to coccolith/biofilm matrix, though the tRNA-synthetase annotation is ambiguous.'),
    'OG0009301': ('Proline-rich-membrane-anchor / SHOCT domain; blast to RNA pol II CTD (likely low-complexity artefact); 11 SP + 11 SP+TM',
                  'weak',
                  'Very high signal-peptide+TM fraction suggests membrane-anchored secreted protein. Function ambiguous; proline-rich repeats resemble RNA-pol CTDs but the localisation pattern is consistent with coccolith-vesicle membrane scaffolds.'),
    'OG0009816': ('α-1,3-mannosyltransferase (MNN13/14/15 family, PF11051); 3 SP / 6 TM',
                  'moderate',
                  'Golgi mannosyltransferase — glycosylates coccolith/matrix glycoproteins. Relevant because biomineralisation matrix proteins are heavily glycosylated and N/O-glycan patterns influence CaCO3 nucleation.'),
    'OG0010177': ('Ankyrin-repeat protein (24 ANK/KRIT1 domains); globular',
                  'unlikely',
                  'Intracellular scaffolding; no direct calcification link.'),
    'OG0010867': ('Pentapeptide-repeat protein (27 PF00805 hits)',
                  'moderate',
                  'See OG0001976 — pentapeptide repeats can carry Ca2+-binding/CCM function in cyanobacteria/algae; plausible indirect role in carbon-concentrating mechanism.'),
    'OG0010991': ('Uncharacterised (no hit); 4 SP / 1 TM',
                  'weak',
                  'Unknown but signal-peptide bearing — cannot rule out secreted matrix protein.'),
    'OG0011061': ('Exostosin-GT47 glycosyltransferase (PF03016); 9 SP+TM, 18 TM — type-II membrane glycosyltransferase',
                  'moderate',
                  'GT47/exostosin transfers sugars to heparan-sulfate / xyloglucan chains. In animals EXT1/2 are essential for skeletal development (heparan-sulfate proteoglycan biosynthesis drives growth-plate calcification); in algae equivalents make extracellular polysaccharide matrix that templates CaCO3. High TM fraction is consistent with Golgi anchoring.'),
    'OG0011197': ('2OG-Fe(II) oxygenase / RAP-domain protein',
                  'weak',
                  '2OG-dioxygenases include collagen prolyl/lysyl hydroxylases (mineralisation-associated in vertebrates), but BLAST points to chloroplast RAP protein. Indirect at best.'),
    'OG0013037': ('Serine aminopeptidase S33 (α/β-hydrolase); blast to titin/obscurin (likely coiled-coil misalignment); 8 SP / 5 TM',
                  'weak',
                  'Likely membrane-anchored serine hydrolase. No direct calcification link; titin hits are generic coiled-coil cross-matches.'),
    'OG0014155': ('Inosine-uridine nucleoside hydrolase (PF01156); 14/22 with signal peptide',
                  'moderate',
                  'High SP fraction is striking for a nucleoside hydrolase — suggests a secreted/periplasmic purine-salvage enzyme possibly active in the coccolith vesicle. Not a core calcification enzyme.'),
    'OG0014250': ('Amidase / Glu-tRNA(Gln) amidotransferase subunit A (PF01425); 17/22 with SP',
                  'weak',
                  'Secretion-signal prevalence hints at extracellular amidase, but no specific calcification role.'),
    'OG0015153': ('Subtelomeric hrmA-cluster protein / piccolo-like; mostly TM or SP+TM',
                  'weak',
                  'Piccolo matches are generic low-complexity; the hrmA-cluster membership suggests a lineage-specific membrane protein — functionally unresolved.'),
    'OG0016203': ('RCC1-repeat protein with a Collagen triple-helix domain; 13/18 TM, 2 SP+TM',
                  'moderate',
                  'Collagen triple-helix modules are mineralisation-associated scaffolds (bone, cartilage, sponge spicules); fusion with RCC1 (Ran-GEF) is unusual. Presence of TM+SP strongly suggests an extracellular collagen-like scaffold.'),
    'OG0016211': ('SUR7/PalI 4-TM membrane family (PF06687); 18/18 TM',
                  'weak',
                  'Fungal SUR7 marks membrane compartments of plasma membrane (MCC/eisosomes). Not classically tied to calcification, but membrane micro-domain organisation could matter for vesicle-based biomineralisation.'),
    'OG0017138': ('Arylsulfatase B (PF00884); 3/16 SP, globular',
                  'moderate',
                  'Arylsulfatases desulfate glycosaminoglycans and polysaccharide matrices. In biomineralising algae, sulfated polysaccharides coat coccoliths; ARSB activity could remodel the organic template. Sulfatases are highly induced during coccolith formation in Emiliania.'),
    'OG0017238': ('tRNA ligase 1',
                  'unlikely',
                  'Housekeeping translation component.'),
    'OG0017305': ('Stf0 sulphotransferase (PF09037); 9/16 SP',
                  'moderate',
                  'Carbohydrate/lipid sulphotransferase. Sulfated polysaccharides and glycolipids are key organic components of coccoliths and algal CaCO3 matrices. High SP fraction consistent with secretion to matrix.'),
    'OG0017361': ('DJ-1/PfpI family cysteine hydrolase / glyoxalase',
                  'unlikely',
                  'Redox/oxidative-stress enzyme; no calcification link.'),
    'OG0017362': ('Red chlorophyll catabolite reductase (PF06405) — blast to AP endonuclease',
                  'unlikely',
                  'Plastid chlorophyll catabolism; no link to calcification.'),
    'OG0017884': ('GNAT acetyltransferase',
                  'unlikely',
                  'Generic acyl-transferase.'),
    'OG0017903': ('FAD-binding oxidoreductase / globin',
                  'unlikely',
                  'O2/redox-handling protein. Not calcification-specific.'),
    'OG0017965': ('Kelch-repeat / β-propeller + NANM (N-acetyl-neuraminate epimerase) — blast to KEAP1/KLHL36 and plant Kelch',
                  'weak',
                  'BTB-Kelch / NANM fusion — substrate-recognition subunit of Cullin E3 ligase with a sialic-acid epimerase module. Sialic acid metabolism can affect mucin/matrix glycosylation but relevance to calcification is indirect.'),
    'OG0018519': ('Uncharacterised; 1 SP / 13 TM',
                  'weak',
                  'Multi-pass membrane protein of unknown function. Possible transporter — cannot rule out ion/HCO3- transport.'),
    'OG0018986': ('FG-GAP / integrin-α-like β-propeller; blast to metabotropic glutamate / sweet-taste 7TM receptor; 8 SP+TM, 2 SP, 2 TM',
                  'strong',
                  'FG-GAP repeats form a Ca2+-binding β-propeller (DxDxDGxxD motifs), the Ca2+-coordinating module of integrin α-subunits. Combined with a 7TM receptor body, this resembles a Ca2+-sensing receptor — directly relevant to Ca2+ uptake regulation during calcification.'),
    'OG0019067': ('Leucine-rich-repeat / TMOD-LMOD-LRR; 3 SP+TM, 2 TM',
                  'weak',
                  'LRR scaffold, possibly membrane-anchored; tropomodulin-like LRR is actin-binding. No direct calcification link.'),
    'OG0019078': ('Glutathione-S-transferase (EF-1γ)',
                  'unlikely',
                  'Cytosolic detoxification enzyme.'),
    'OG0019174': ('Glu-tRNA(Gln) amidotransferase subunit A (PF01425 amidase); 12/13 SP',
                  'weak',
                  'Nearly universal secretion signal is striking — suggests extracellular amidase activity, but no dedicated calcification role.'),
    'OG0019180': ('Uncharacterised',
                  'weak',
                  'Uncharacterised; cannot be assessed.'),
    'OG0019790': ('HECT E3 ubiquitin ligase with ankyrin repeats',
                  'unlikely',
                  'Ubiquitin-mediated protein turnover.'),
    'OG0019817': ('FAD-binding monooxygenase / zeaxanthin epoxidase',
                  'unlikely',
                  'Carotenoid / redox metabolism.'),
    'OG0019873': ('Class-I/II aminotransferase (1-aminocyclopropane-1-carboxylate synthase)',
                  'unlikely',
                  'Amino-acid/ethylene metabolism.'),
    'OG0019887': ('Mitochondrial ribosome-associated GTPase (MMR_HSR1) / intron maturase',
                  'unlikely',
                  'Mitoribosome biogenesis.'),
    'OG0019917': ('Cyclic-nucleotide-binding / Apolipoprotein A domain; 7 SP, 2 SP+TM',
                  'weak',
                  'cNMP-binding regulatory subunit; high SP unusual — possibly a secreted ligand-binding module. No direct calcification evidence.'),
    'OG0020622': ('AAA+ RecD-like DNA helicase',
                  'unlikely',
                  'DNA repair.'),
    'OG0020657': ('Uncharacterised; 7/11 SP',
                  'weak',
                  'Putative secreted protein, function unknown — could be a novel lineage-specific matrix protein.'),
    'OG0020696': ('Uncharacterised; 3 TM, 1 SP',
                  'weak',
                  'Unknown multi-pass membrane protein.'),
    'OG0020700': ('PDZ-domain scaffold',
                  'unlikely',
                  'Scaffolding; no known calcification tie-in.'),
    'OG0020703': ('Heparan-sulfate 6-O-sulfotransferase (PF13469 / PF00685); 6/11 SP',
                  'strong',
                  '6-O-sulfation of heparan-sulfate and related glycans is essential for assembly of biomineralisation matrix (in animals HS6ST is required for skeletal patterning; in coccolithophores sulfated polysaccharides form the coccolith base plate).'),
    'OG0020706': ('Uncharacterised; 4 TM, 1 SP, 1 SP+TM',
                  'weak',
                  'Multi-pass membrane protein of unknown function; transport function cannot be excluded.'),
    'OG0020728': ('Uncharacterised',
                  'weak', 'No data.'),
    'OG0020729': ('Uncharacterised',
                  'weak', 'No data.'),
    'OG0021347': ('Trypsin-family serine protease (PF00089); 10/10 SP — fully secreted',
                  'moderate',
                  '100 % signal-peptide — clearly secreted protease. Extracellular proteases shape the organic matrix around coccoliths/spicules; involvement in matrix remodelling during calcification is plausible.'),
    'OG0021406': ('Short-chain dehydrogenase/reductase',
                  'unlikely', 'Generic SDR metabolic enzyme.'),
    'OG0021517': ('F-box protein (SCF ubiquitin-ligase adaptor)',
                  'unlikely', 'Ubiquitin-dependent protein degradation.'),
    'OG0021520': ('Uncharacterised; 1 BETA + 9 GLOB',
                  'weak',
                  'One predicted β-barrel outer-membrane-like sequence — lineage-specific.'),
    'OG0021523': ('Protein phosphatase 2C (PP2C)',
                  'weak',
                  'Ser/Thr phosphatase — Ca2+-independent; only tangential calcification link via signalling.'),
    'OG0021531': ('PDZ + GRASP55/65 + S41 peptidase; 1 SP',
                  'weak',
                  'GRASP55/65 tethers Golgi membranes; possible role in Golgi-derived coccolith-vesicle trafficking.'),
    'OG0022355': ('DUF2237 bacterial-like protein; 6 SP+TM, 2 SP',
                  'weak', 'Unknown; secreted/membrane but uncharacterised.'),
    'OG0022381': ('FAD/NAD(P)H oxidoreductase / flavodoxin',
                  'unlikely', 'Redox enzyme.'),
    'OG0022455': ('Uncharacterised; 2 SP', 'weak', 'Unknown.'),
    'OG0022473': ('Uncharacterised; 3 SP+TM, 1 SP', 'weak', 'Secreted/membrane-anchored unknown protein — possible novel matrix component.'),
    'OG0022474': ('F-box / WD-repeat adaptor', 'unlikely', 'Ubiquitination.'),
    'OG0022492': ('β-propeller / WD40 (TEP1/WDR3/WDR5/WDR55/THOC3-like)',
                  'unlikely', 'Ribosome/spliceosome biogenesis WD40 scaffolds.'),
    'OG0022500': ('Uncharacterised', 'weak', 'Unknown.'),
    'OG0022520': ('Uncharacterised', 'weak', 'Unknown.'),
    'OG0022524': ('Pentapeptide-repeat protein',
                  'moderate', 'See OG0001976 — potential CCM / Ca2+-binding cyanobacterial-type repeat.'),
    'OG0022528': ('Uncharacterised', 'weak', 'Unknown.'),
    'OG0023496': ('Protein kinase (blast: BMP-2-inducible kinase / AP2-associated kinase)',
                  'moderate',
                  'BMP2K/AAK1 homologues; BMP-2 signalling is central to osteoblast differentiation and vertebrate bone mineralisation. Kinase activity could be involved in analogous calcification signalling.'),
    'OG0023566': ('RHD3/SEY1 ER-membrane-fusion GTPase',
                  'moderate',
                  'SEY1/RHD3 drives ER-network homotypic fusion. ER remodelling is upstream of coccolith-vesicle formation; disruption of SEY1 alters secretion — relevant to biomineralisation machinery.'),
    'OG0023594': ('Glycosyl-transferase family 8 (PF01501); 1 SP+TM',
                  'moderate',
                  'GT8 includes galacturonosyl/galactosyltransferases that build pectic and other cell-wall/matrix polysaccharides; in algae they contribute to the organic coccolith matrix.'),
    'OG0023657': ('Uncharacterised', 'weak', 'Unknown.'),
    'OG0023790': ('Red chlorophyll catabolite reductase',
                  'unlikely', 'Chlorophyll catabolism.'),
    'OG0024846': ('TOD1/MUCI70 glycosyltransferase-like domain (PF04765); 2 SP',
                  'moderate',
                  'MUCI70 is a mucin-type O-glycosyltransferase required for cell-wall biosynthesis in plants; in algae the homologue likely glycosylates matrix mucin-like proteins forming the coccolith organic scaffold.'),
    'OG0024979': ('Uncharacterised', 'weak', 'Unknown.'),
    'OG0025015': ('Uncharacterised; 7/7 TM',
                  'weak',
                  'All predicted as multi-pass membrane protein — possible novel transporter; function unassigned.'),
    'OG0025016': ('Uncharacterised', 'weak', 'Unknown.'),
    'OG0025059': ('Uncharacterised; 5/7 SP',
                  'weak',
                  'Secretion-prone unknown protein — possible lineage-specific matrix protein.'),
    'OG0026865': ('JmjC-domain Fe(II)/2OG hydroxylase',
                  'weak',
                  'JmjC hydroxylases include histone demethylases and the PSR arginine-demethylase. A subset (2OG-dioxygenases) performs collagen Pro/Lys hydroxylation, coupling to mineralisation, but BLAST here points to demethylase activity.'),
}

# ---------------------------------------------------------------------------
# 8. Build output rows
# ---------------------------------------------------------------------------
rows = []
for og in ogs:
    n_seq = len(og_seqs[og])
    tm = og_tm[og]
    sp = og_sp[og]
    n_sp = sp.get('Signal Peptide (Sec/SPI)', 0)
    n_tm = tm.get('TM', 0) + tm.get('SP+TM', 0)
    n_glob = tm.get('GLOB', 0)
    tm_breakdown = ', '.join(f'{k}={v}' for k, v in sorted(tm.items()))
    top_b = consensus_blast(og_blast.get(og, []), k=3)
    blast_summary = '; '.join(short_desc(h['desc']) for h in top_b) if top_b else 'no blast hit'
    blast_best_ev = min((h['evalue'] for h in og_blast.get(og, [])), default=None)
    top_p = consensus_pfam(og_pfam.get(og, []), k=6)
    pfam_summary = ('; '.join(f'{p} ({r["acc"]}, n={r["n"]}, {r["desc"]})'
                              for p, r in top_p)
                    if top_p else 'no PFAM hit')
    loc = []
    if n_sp: loc.append(f'{n_sp}/{n_seq} w/ signal peptide')
    if n_tm: loc.append(f'{n_tm}/{n_seq} with TM')
    if n_glob and not loc: loc.append(f'{n_glob}/{n_seq} globular')
    func, relev, note = CURATED.get(og, ('uncurated', 'weak', ''))
    rows.append({
        'OG': og, 'n_seq': n_seq,
        'function_summary': func,
        'top_pfam_domains': pfam_summary,
        'top_blast_hits': blast_summary,
        'blast_best_evalue': f'{blast_best_ev:.1e}' if blast_best_ev is not None else 'NA',
        'signalp_sp_count': n_sp,
        'tmhmm_breakdown': tm_breakdown,
        'localization_note': '; '.join(loc) if loc else 'NA',
        'calcification_relevance': relev,
        'calcification_rationale': note,
    })

# ---------------------------------------------------------------------------
# 9. Write TSV
# ---------------------------------------------------------------------------
out_path = os.path.join(ROOT, 'orthogroup_annotation.tsv')
fields = ['OG', 'n_seq', 'function_summary',
          'top_pfam_domains', 'top_blast_hits', 'blast_best_evalue',
          'signalp_sp_count', 'tmhmm_breakdown', 'localization_note',
          'calcification_relevance', 'calcification_rationale']
with open(out_path, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader()
    for r in rows:
        w.writerow(r)

print(f'Wrote {out_path}')
cc = Counter(r['calcification_relevance'] for r in rows)
print('Calcification relevance counts:', dict(cc))
print('\nStrong / moderate candidates:')
for r in rows:
    if r['calcification_relevance'] in ('strong', 'moderate'):
        print(f"  {r['OG']} [{r['calcification_relevance']}] — {r['function_summary']}")
