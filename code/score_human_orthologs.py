#!/usr/bin/env python3
"""Map rat cancer missense variants to human orthologs and score with CADD/AlphaMissense/SIFT."""
import json, time, urllib.request
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

def get_json(url):
    for a in range(5):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "biomni"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.load(r)
        except Exception as e:
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

def human_canonical_tx(gene):
    d = get_json(f"https://rest.ensembl.org/lookup/symbol/human/{gene}?expand=1")
    if not d:
        return None
    for t in d.get("Transcript", []):
        if t.get("is_canonical"):
            return t["id"]
    return None

def parse_fasta(fa):
    if not fa:
        return None
    lines = fa.splitlines()
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    return "".join(lines)

def ensembl_protein(tx):
    return parse_fasta(get_text(f"https://rest.ensembl.org/sequence/id/{tx}?type=protein;content-type=text/plain"))

def ensembl_cds(tx):
    return parse_fasta(get_text(f"https://rest.ensembl.org/sequence/id/{tx}?type=cds;content-type=text/plain"))

def rat_protein(nm):
    fa = get_text(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={nm}&rettype=fasta_cds_aa&retmode=fasta")
    return parse_fasta(fa)

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = "global"
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

def map_rat_to_human(rat_seq, hum_seq, rat_pos):
    aln = aligner.align(rat_seq, hum_seq)[0]
    tb, qb = aln.aligned  # target=rat, query=human; each (n,2) [start,end)
    for i in range(len(tb)):
        rs, re_ = int(tb[i][0]), int(tb[i][1])
        hs = int(qb[i][0])
        if rs <= rat_pos - 1 < re_:
            return hs + (rat_pos - 1 - rs) + 1
    return None

def resolve_hgvs(cds, pos, ref_aa, alt_aa):
    import itertools
    cs = (pos - 1) * 3
    codon = cds[cs:cs + 3]
    if str(Seq(codon).translate()) != ref_aa:
        return None, codon
    # single-nucleotide path
    for i in range(3):
        for nt in "ACGT":
            if nt == codon[i]:
                continue
            nc = codon[:i] + nt + codon[i + 1:]
            if str(Seq(nc).translate()) == alt_aa:
                return f"c.{cs + i + 1}{codon[i]}>{nt}", codon
    # MNV fallback: whole-codon delins to the alt codon with fewest nt differences
    alt_codons = ["".join(c) for c in itertools.product("ACGT", repeat=3)
                  if str(Seq("".join(c)).translate()) == alt_aa]
    alt_codons.sort(key=lambda c: sum(x != y for x, y in zip(c, codon)))
    best = alt_codons[0]
    diffs = [i for i in range(3) if codon[i] != best[i]]
    s, e = min(diffs), max(diffs)
    g_s, g_e = cs + s + 1, cs + e + 1
    if s == e:
        return f"c.{g_s}{codon[s]}>{best[s]}", codon
    return f"c.{g_s}_{g_e}delins{best[s:e+1]}", codon

AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G",
       "His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S",
       "Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*"}

def parse_hgvsp(h):
    # p.Gly566Asp -> (G, 566, D)
    b = h.replace("p.", "")
    ref = AA3[b[:3]]
    import re
    m = re.search(r"\d+", b)
    pos = int(m.group())
    alt = AA3[b[m.end():m.end()+3]]
    return ref, pos, alt

# (gene, rat_tx, hgvs_p)
TODO = [
    ("AFF3",  "XM_063267369.1", "p.Gly566Asp"),
    ("ARID1B","NM_001419802.1", "p.Ala1037Thr"),
    ("ARID1B","NM_001419802.1", "p.Gly1948Asp"),
    ("CLTC",  "XM_063269683.1", "p.Ser1650Asn"),
    ("KMT2D", "XM_063262838.1", "p.Gly1933Ser"),
    ("KMT2D", "XM_063262838.1", "p.Gly1866Ser"),
    ("TET1",  "XM_039099324.2", "p.Val926Ile"),
    ("TP53",  "NM_030989.3",    "p.Arg271Ser"),
]

# cache proteins per gene
cache = {}
rows = []
for gene, rat_tx, hgvsp in TODO:
    ref_aa, rat_pos, alt_aa = parse_hgvsp(hgvsp)
    if gene not in cache:
        htx = human_canonical_tx(gene)
        hprot = ensembl_protein(htx)
        hcds = ensembl_cds(htx)
        rprot = rat_protein(rat_tx)
        cache[gene] = (htx, hprot, hcds, rprot)
        time.sleep(0.3)
    htx, hprot, hcds, rprot = cache[gene]
    if not all([htx, hprot, hcds, rprot]):
        rows.append((gene, hgvsp, None, None, "FETCH_FAIL", None, None, None, None, None))
        continue
    # verify rat ref aa
    rat_ref_ok = rprot[rat_pos-1] == ref_aa
    hpos = map_rat_to_human(rprot, hprot, rat_pos)
    if hpos is None:
        rows.append((gene, hgvsp, None, None, "NOT_CONSERVED", None, None, None, None, None))
        continue
    hum_ref = hprot[hpos-1]
    if hum_ref != ref_aa:
        rows.append((gene, hgvsp, hpos, hum_ref, "REF_MISMATCH", None, None, None, None, None))
        continue
    hgvs_c, codon = resolve_hgvs(hcds, hpos, ref_aa, alt_aa)
    if hgvs_c is None:
        rows.append((gene, hgvsp, hpos, hum_ref, "CODON_FAIL", None, None, None, None, None))
        continue
    # human VEP
    v = get_json(f"https://rest.ensembl.org/vep/human/hgvs/{htx}:{hgvs_c}?content-type=application/json;CADD=1;AlphaMissense=1;SIFT=1")
    cadd = am_class = am_score = hsift = None
    if v:
        for tc in v[0].get("transcript_consequences", []):
            if tc.get("transcript_id") == htx:
                cadd = tc.get("cadd_phred")
                hsift = tc.get("sift_prediction")
                am = tc.get("alphamissense") or {}
                am_class = am.get("am_class")
                am_score = am.get("am_pathogenicity")
    hum_hgvsp = f"p.{ref_aa}{hpos}{alt_aa}"
    rows.append((gene, hgvsp, hpos, hum_ref, "OK", hum_hgvsp, cadd, am_class, am_score, hsift))
    print(f"{gene} {hgvsp} -> human {hum_hgvsp} ({htx} {hgvs_c}): CADD={cadd} AM={am_class}({am_score}) SIFT={hsift}")
    time.sleep(0.3)

import pandas as pd
df = pd.DataFrame(rows, columns=["gene","rat_hgvsp","human_pos","human_ref","status",
                                  "human_hgvsp","cadd_phred","am_class","am_score","human_sift"])
df.to_csv("/workspace/hede_followup/variants/human_ortholog_scores.tsv", sep="\t", index=False)
print("\nwrote variants/human_ortholog_scores.tsv")
print(df.to_string())
