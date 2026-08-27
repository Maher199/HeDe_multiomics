#!/usr/bin/env python3
"""Confirm Cdkn2a deletion breakpoints at base resolution from split reads.
Sniffles2 DEL.1A18S5: NC_086023.1:109021504 -> 109151061 (SVLEN=-129557).
For each split read (SA tag), recover the pair of alignment blocks, find the
read-relative junction, and map it to reference coordinates on both sides.
Microhomology = reference bases shared at the two breakpoint ends.
"""
import pysam, re
from collections import Counter

BAM = "/mnt/shared-workspace/shared/HeDe_on_GRCr8.sorted.bam"
FA  = "/mnt/shared-workspace/shared/rn_GRCr8.fa"
CHROM = "NC_086023.1"
# Sniffles2 call
BP_L, BP_R = 109021504, 109151061   # 1-based POS / END
REG = (CHROM, 109014000, 109158500)

bam = pysam.AlignmentFile(BAM, "rb")
fa  = pysam.FastaFile(FA)

def parse_sa(sa):
    out = []
    for e in sa.rstrip(";").split(";"):
        p = e.split(",")
        if len(p) >= 6:
            out.append(dict(rname=p[0], pos=int(p[1]), strand=p[2],
                            cigar=p[3], mapq=int(p[4]), nm=int(p[5])))
    return out

def cig_ref_consumed(cigar):
    return sum(l for op, l in cigar if op in (0, 2, 3, 7, 8))  # M D N = X

def cig_query_left_clip(cigar):
    return cigar[0][1] if cigar and cigar[0][0] in (4, 5) else 0

def cig_query_right_clip(cigar):
    return cigar[-1][1] if cigar and cigar[-1][0] in (4, 5) else 0

# Collect split-read junctions that bridge the two breakpoints
junct = []   # (left_ref_end_1based, right_ref_start_1based, readname, nSA)
for r in bam.fetch(*REG):
    if r.is_secondary or r.is_supplementary:
        continue
    if not r.has_tag("SA"):
        continue
    # primary block
    prim = dict(rname=r.reference_name, pos=r.reference_start + 1,  # 1-based
                strand="-" if r.is_reverse else "+",
                cigar=r.cigarstring, mapq=r.mapping_quality,
                ref_start=r.reference_start + 1,
                ref_end=r.reference_end,  # 1-based inclusive-ish
                lc=cig_query_left_clip(r.cigartuples),
                rc=cig_query_right_clip(r.cigartuples),
                qlen=r.query_length)
    blocks = [prim] + [dict(b, ref_start=b["pos"],
                            ref_end=None, lc=None, rc=None, qlen=r.query_length)
                       for b in parse_sa(r.get_tag("SA"))]
    # compute ref_end and clips for SA blocks from their cigar
    for b in blocks[1:]:
        ct = pysam.AlignedSegment().cigarstring  # placeholder
        m = re.findall(r"(\d+)([MIDNSHP=X])", b["cigar"])
        ops = [(int(l), o) for l, o in m]
        b["ref_end"] = b["pos"] + sum(l for l, o in ops if o in "MDN=X") - 1
        b["lc"] = ops[0][1] if ops and ops[0][1] in "SH" else 0
        b["rc"] = ops[-1][1] if ops and ops[-1][1] in "SH" else 0
    # does any block end near BP_L and another start near BP_R?
    ends_L = [b for b in blocks if abs(b["ref_end"] - BP_L) <= 50]
    starts_R = [b for b in blocks if abs(b["pos"] - (BP_R + 1)) <= 50]
    if ends_L and starts_R:
        for bl in ends_L:
            for br in starts_R:
                junct.append((bl["ref_end"], br["pos"], r.query_name, len(blocks) - 1))

print(f"split reads bridging the two breakpoints: {len(junct)}")
le = Counter(j[0] for j in junct)
rs = Counter(j[1] for j in junct)
print("\nLeft-breakpoint (ref end of left block), top values:")
for pos, n in le.most_common(6):
    print(f"  {pos}  ({n} reads)   [Sniffles2 POS={BP_L}]")
print("\nRight-breakpoint (ref start of right block), top values:")
for pos, n in rs.most_common(6):
    print(f"  {pos}  ({n} reads)   [Sniffles2 END+1={BP_R+1}]")

# Consensus breakpoints
if le and rs:
    cl = le.most_common(1)[0][0]   # last reference base kept on left
    cr = rs.most_common(1)[0][0]   # first reference base kept on right
    print(f"\nConsensus: left block ends at {cl}, right block starts at {cr}")
    print(f"Deleted interval (reference): {cl+1} .. {cr-1}  (len {cr-1-(cl+1)+1} bp)")
    print(f"Sniffles2:                    {BP_L+1} .. {BP_R}  (len {BP_R-BP_L} bp)")

    # Microhomology: compare reference 3' of left junction vs 5' of right junction.
    # The junction joins left-cl .. right-cr. Microhomology = longest suffix of the
    # left-kept flank that equals a prefix of the right-kept flank (ambiguous assignment).
    L = 30
    left_seq  = fa.fetch(CHROM, cl - L, cl)        # 0-based half-open: last L bases of left flank
    right_seq = fa.fetch(CHROM, cr - 1, cr - 1 + L)  # first L bases of right flank
    left_seq = left_seq.upper(); right_seq = right_seq.upper()
    mh = 0
    for k in range(min(len(left_seq), len(right_seq)), 0, -1):
        if left_seq[-k:] == right_seq[:k]:
            mh = k; break
    print(f"\nMicrohomology at junction: {mh} bp")
    print(f"  left flank (..{cl}):  ...{left_seq}")
    print(f"  right flank ({cr}..): ...{right_seq}")
    if mh > 0:
        print(f"  shared {mh} bp: {left_seq[-mh:]}")
    mech = ("MMEJ (microhomology-mediated)" if mh >= 2 else
            "NHEJ (blunt/1-bp)" if mh <= 1 else "NHEJ")
    print(f"  inferred repair mechanism: {mech}")

bam.close(); fa.close()
print("\nDONE")
