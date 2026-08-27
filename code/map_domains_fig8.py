#!/usr/bin/env python3
"""Fetch human UniProt domains, map to rat protein, then to GRCr8 genomic coords (fig8 second-tier drivers)."""
import json, time, urllib.request, os
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

WD = "/workspace/hede_followup"
OUT = os.path.join(WD, "fig8_data")

GENES = {  # gene -> (chrom, strand, rat_tx, human_uniprot)
    "KMT2C": ("NC_086022.1", "+", "XM_039108963.2", "Q8NEZ4"),
    "ARID1B":("NC_086019.1", "+", "NM_001419802.1", "Q8NFD5"),
    "TET1":  ("NC_086038.1", "+", "XM_039099324.2", "Q8NFU7"),
    "CLTC":  ("NC_086028.1", "-", "XM_063269683.1", "Q00610"),
    "WNK2":  ("NC_086035.1", "+", "XM_006253715.5", "Q9Y3S1"),
    "DHX9":  ("NC_086031.1", "-", "NM_001395556.1", "Q08211"),
}
FEAT_TYPES = {"Domain", "Region", "Repeat", "Zinc finger", "DNA binding", "Motif", "Active site", "Binding site"}

def get_json(url):
    for a in range(5):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "biomni"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.load(r)
        except Exception:
            time.sleep(1.5 * (a + 1))
    return None

def get_text(url):
    for a in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "biomni"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read().decode()
        except Exception:
            time.sleep(1.5 * (a + 1))
    return None

def parse_fasta(fa):
    lines = fa.splitlines()
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    return "".join(lines)

def rat_protein(nm):
    return parse_fasta(get_text(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={nm}&rettype=fasta_cds_aa&retmode=fasta"))

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = "global"
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

def build_h2r_map(rat_seq, hum_seq):
    aln = aligner.align(rat_seq, hum_seq)[0]
    tb, qb = aln.aligned
    m = {}
    for i in range(len(tb)):
        rs, re_ = int(tb[i][0]), int(tb[i][1])
        hs = int(qb[i][0])
        for k in range(re_ - rs):
            m[hs + k + 1] = rs + k + 1
    return m

def aa_to_genomic_map(cds_df, strand):
    exons = sorted([(int(r.start), int(r.end)) for r in cds_df.itertuples()])
    nucs = []
    if strand == "+":
        for s, e in exons:
            nucs.extend(range(s, e + 1))
    else:
        for s, e in sorted(exons, reverse=True):
            nucs.extend(range(e, s - 1, -1))
    return nucs

def aa_range_to_genomic(nucs, a, b):
    i0 = max(0, (a - 1) * 3)
    i1 = min(len(nucs) - 1, (b - 1) * 3 + 2)
    if i0 >= len(nucs):
        return None
    g = nucs[i0:i1 + 1]
    return min(g), max(g)

for gene, (chrom, strand, rat_tx, up) in GENES.items():
    cds_df = pd.read_csv(f"{OUT}/cds_{gene}.tsv", sep="\t")
    nucs = aa_to_genomic_map(cds_df, strand)
    rprot = rat_protein(rat_tx)
    u = get_json(f"https://rest.uniprot.org/uniprotkb/{up}.json")
    hprot = u["sequence"]["value"]
    h2r = build_h2r_map(rprot, hprot)
    rows = []
    for ft in u.get("features", []):
        if ft.get("type") not in FEAT_TYPES:
            continue
        desc = ft.get("description", ft.get("type"))
        a = int(ft["location"]["start"]["value"])
        b = int(ft["location"]["end"]["value"])
        ra = h2r.get(a)
        rb = h2r.get(b)
        if ra is None or rb is None:
            mapped = [h2r[p] for p in range(a, b + 1) if p in h2r]
            if not mapped:
                continue
            ra, rb = min(mapped), max(mapped)
        gr = aa_range_to_genomic(nucs, min(ra, rb), max(ra, rb))
        if gr is None:
            continue
        rows.append((gene, ft.get("type"), desc, a, b, min(ra, rb), max(ra, rb), gr[0], gr[1]))
    df = pd.DataFrame(rows, columns=["gene", "ftype", "description", "hum_start", "hum_end",
                                     "rat_start", "rat_end", "gstart", "gend"])
    df.to_csv(f"{OUT}/domains_{gene}.tsv", sep="\t", index=False)
    print(f"{gene} ({up}): {len(df)} features; rat_prot={len(rprot)}aa hum_prot={len(hprot)}aa")
    for _, r in df.iterrows():
        print(f"    {r.ftype:12s} {str(r.description)[:34]:34s} hum {r.hum_start}-{r.hum_end} -> rat {r.rat_start}-{r.rat_end} -> g {r.gstart}-{r.gend}")
    time.sleep(0.3)
print("done")
