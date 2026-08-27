#!/usr/bin/env python3
"""Read-level phasing of Tp53 R273S (C>A) with nearby heterozygous sites.
If two haplotypes are present (one carrying A at R273S, one carrying C), the
locus is heterozygous with a retained WT allele (no LOH). Long HiFi reads let a
single read span R273S and a neighbouring het site, giving direct cis-phasing.
All coordinates converted to 0-based for pysam.
"""
import pysam
from collections import Counter, defaultdict

BAM = "/mnt/shared-workspace/shared/HeDe_on_GRCr8.sorted.bam"
CHROM = "NC_086028.1"
R273S_0 = 54807960          # VCF POS 54807961, C>A
R_REF, R_ALT = "C", "A"
# nearby het sites (VCF POS 1-based, REF, ALT); all simple deletions (REF len>ALT len)
SITES = [
    ("i54802512", 54802512, "Ca", "C"),
    ("i54808916", 54808916, "Ga", "G"),
    ("i54811961", 54811961, "Ga", "G"),
    ("i54812994", 54812994, "Ga", "G"),
    ("i54813761", 54813761, "Ca", "C"),
]

bam = pysam.AlignmentFile(BAM, "rb")

def refmap(read):
    m = {}
    for q, r in read.get_aligned_pairs():
        if r is not None:
            m[r] = q
    return m

def snv_allele(read, m, pos0, ref, alt):
    q = m.get(pos0)
    if q is None:
        return None
    b = read.query_sequence[q].upper()
    if b == ref: return "REF"
    if b == alt: return "ALT"
    return None

def del_allele(read, m, pos0):
    # anchor base at pos0 (VCF POS-1); deleted base at pos0+1
    if pos0 not in m or (pos0 + 1) not in m:
        return None
    return "REF" if m[pos0 + 1] is not None else "ALT"

r273s_counts = Counter()
pair_tab = defaultdict(Counter)   # site -> (r273s_allele, site_allele) -> n
n_reads = 0
for read in bam.fetch(CHROM, R273S_0 - 9000, R273S_0 + 9000):
    if read.is_secondary or read.is_supplementary or read.mapping_quality < 20:
        continue
    m = refmap(read)
    ra = snv_allele(read, m, R273S_0, R_REF, R_ALT)
    if ra is None:
        continue
    n_reads += 1
    r273s_counts[ra] += 1
    for name, pos1, ref, alt in SITES:
        sa = del_allele(read, m, pos1 - 1)
        if sa is not None:
            pair_tab[name][(ra, sa)] += 1
bam.close()

print(f"reads spanning R273S (MAPQ>=20): {n_reads}")
print(f"R273S read-level alleles: REF(C)={r273s_counts['REF']}  ALT(A)={r273s_counts['ALT']}"
      f"  -> alt fraction {r273s_counts['ALT']/max(1,n_reads):.3f}")
print("\nPairwise cis-phasing (R273S allele x nearby-site allele):")
for name, pos1, ref, alt in SITES:
    t = pair_tab[name]
    tot = sum(t.values())
    if tot == 0:
        continue
    print(f"\n  {name} (POS {pos1}, {ref}>{alt}), {tot} informative reads:")
    for k in [("REF","REF"),("REF","ALT"),("ALT","REF"),("ALT","ALT")]:
        print(f"    R273S={k[0]:3} site={k[1]:3} : {t.get(k,0)}")
    # haplotype inference
    rr, ra_, ar, aa = t.get(("REF","REF"),0), t.get(("REF","ALT"),0), t.get(("ALT","REF"),0), t.get(("ALT","ALT"),0)
    # linkage: does R273S-ALT travel with site-REF or site-ALT?
    if (ar + aa) > 0 and (rr + ra_) > 0:
        alt_with = "site-REF" if ar > aa else ("site-ALT" if aa > ar else "unclear")
        ref_with = "site-REF" if rr > ra_ else ("site-ALT" if ra_ > rr else "unclear")
        print(f"    -> R273S-ALT(A) travels with {alt_with}; R273S-REF(C) travels with {ref_with}")
print("\nDONE")
