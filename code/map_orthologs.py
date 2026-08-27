#!/usr/bin/env python3
"""Map rat RefSeq protein positions to human ortholog positions via global alignment."""
from Bio import Align
from Bio.Align import substitution_matrices

def read_fasta(path):
    seq = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq

pairs = {
    "DICER1": ("protein_seqs/rat_NM_001427215.1.fasta", "protein_seqs/human_Q9UPY3.fasta", [1701]),
    "FAT1":   ("protein_seqs/rat_XM_063275777.1.fasta", "protein_seqs/human_Q14517.fasta", [2257]),
    "POLQ":   ("protein_seqs/rat_NM_001105878.1.fasta", "protein_seqs/human_O75417.fasta", [1509, 1510]),
}

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

for gene, (rat_path, hum_path, positions) in pairs.items():
    rat = read_fasta(rat_path)
    hum = read_fasta(hum_path)
    aln = aligner.align(hum, rat)[0]   # target=human, query=rat
    # aln.aligned -> (blocks_target(human), blocks_query(rat)); each tuple of (start,end)
    hblocks, rblocks = aln.aligned
    # build rat(1-based) -> human(1-based) map
    r2h = {}
    for (hs, he), (rs, re_) in zip(hblocks, rblocks):
        for k in range(re_ - rs):
            r2h[rs + k + 1] = hs + k + 1   # convert to 1-based
    ident = sum(1 for (hs,he),(rs,re_) in zip(hblocks,rblocks) for k in range(re_-rs) if hum[hs+k]==rat[rs+k])
    alnlen = sum(he-hs for hs,he in hblocks)
    print(f"\n=== {gene} ===  rat_len={len(rat)} human_len={len(hum)} aligned_cols={alnlen} identity={ident/alnlen*100:.1f}%")
    for pos in positions:
        hpos = r2h.get(pos)
        raa = rat[pos-1]
        if hpos is None:
            print(f"  rat {raa}{pos} -> human: GAP (no aligned residue)")
        else:
            haa = hum[hpos-1]
            match = "IDENTICAL" if haa==raa else f"DIFFERS (human ref={haa})"
            print(f"  rat {raa}{pos} -> human {haa}{hpos}  [{match}]")
