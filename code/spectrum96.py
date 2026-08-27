#!/usr/bin/env python3
"""Extract the 96-channel trinucleotide spectrum from a somatic SNV VCF.

Mirrors the PacBio analysis: channel format "A[C>A]A", pyrimidine-reference
orientation (C>T/G>A etc. folded), 100% centre-base reference-match validation
(hard-fails if any SNV's REF disagrees with the FASTA).

Usage: spectrum96.py <vcf(.gz)> <ref.fa> <out_prefix>
Writes <out_prefix>.96spectrum.csv (channel,count) and <out_prefix>.6class.csv.
"""
import sys
import gzip
import csv
from collections import Counter

vcf_path, ref_path, out_prefix = sys.argv[1], sys.argv[2], sys.argv[3]

seqs = {}
name = None
buf = []
with open(ref_path) as f:
    for line in f:
        if line.startswith(">"):
            if name is not None and name.startswith("NC_"):
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
if name is not None and name.startswith("NC_"):
    seqs[name] = "".join(buf)
print(f"loaded {len(seqs)} chromosome sequences", file=sys.stderr)

COMP = str.maketrans("ACGT", "TGCA")
def rc(s):
    return s.translate(COMP)[::-1]

BASES = ["A", "C", "G", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
channels = [f"{l}[{s}]{r}" for s in SUBS for l in BASES for r in BASES]

spec = Counter()
six = Counter()
n_snv = n_pass = n_refmatch = 0
skipped = Counter()

opener = gzip.open if vcf_path.endswith(".gz") else open
with opener(vcf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        n_snv += 1
        p = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt, filt = p[0], int(p[1]), p[3].upper(), p[4].upper(), p[6]
        if filt not in ("PASS", "."):
            skipped["nonPASS"] += 1
            continue
        alts = alt.split(",")
        if len(alts) != 1 or len(ref) != 1 or len(alts[0]) != 1:
            skipped["not_biallelic_snv"] += 1
            continue
        if ref not in "ACGT" or alts[0] not in "ACGT":
            skipped["non_acgt"] += 1
            continue
        seq = seqs.get(chrom)
        if seq is None:
            skipped["scaffold"] += 1
            continue
        i = pos - 1
        if i <= 0 or i >= len(seq) - 1:
            skipped["edge"] += 1
            continue
        tri = seq[i - 1:i + 2].upper()
        if tri[1] != ref:
            print(f"FATAL: REF mismatch at {chrom}:{pos} VCF={ref} FASTA={tri[1]}", file=sys.stderr)
            sys.exit(1)
        n_refmatch += 1
        n_pass += 1
        sub = f"{ref}>{alts[0]}"
        if ref in "CT":
            channel = f"{tri[0]}[{sub}]{tri[2]}"
        else:
            channel = f"{rc(tri)[0]}[{rc(ref)}>{rc(alts[0])}]{rc(tri)[2]}"
            sub = f"{rc(ref)}>{rc(alts[0])}"
        spec[channel] += 1
        six[sub] += 1

print(f"VCF records: {n_snv}; used SNVs: {n_pass}; REF==FASTA: {n_refmatch} (100% required)",
      file=sys.stderr)
print(f"skipped: {dict(skipped)}", file=sys.stderr)

with open(out_prefix + ".96spectrum.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["channel", "count"])
    for ch in channels:
        w.writerow([ch, spec.get(ch, 0)])
with open(out_prefix + ".6class.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["subtype", "count"])
    for s in SUBS:
        w.writerow([s, six.get(s, 0)])
print(f"wrote {out_prefix}.96spectrum.csv and .6class.csv", file=sys.stderr)
