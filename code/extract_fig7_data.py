#!/usr/bin/env python3
"""Extract per-gene data for fig7 browser views: coverage, gene model, mutations."""
import subprocess, gzip, re, os
import pandas as pd

WD = "/workspace/hede_followup"
BAM = "/mnt/shared-workspace/shared/HeDe_on_GRCr8.sorted.bam"
GTF = "/workspace/hede_followup/ref/GCF_036323735.1_GRCr8_genomic.gtf.gz"
SAM = "/workspace/.conda/envs/persist/bin/samtools"
OUT = os.path.join(WD, "fig7_data")
os.makedirs(OUT, exist_ok=True)

# gene -> (chrom, start, end, strand, transcript)
GENES = {
    "TP53":  ("NC_086028.1", 54798871, 54810300, "+", "NM_030989.3"),
    "KMT2D": ("NC_086025.1", 131859696, 131901032, "-", "XM_063262838.1"),
    "EP300": ("NC_086025.1", 114987857, 115058652, "+", "XM_039080287.2"),
    "FAT1":  ("NC_086034.1", 53909759, 54029175, "-", "XM_063275777.1"),
    "DICER1":("NC_086024.1", 129392298, 129457252, "-", "NM_001427215.1"),
}

# --- mutations ---
cc = pd.read_csv(f"{WD}/variants/cancer_coding_mutations.tsv", sep="\t")
cc = cc.drop_duplicates(subset=["human_symbol", "CHROM", "POS"])
mut = cc[cc["human_symbol"].isin(GENES)].copy()
mut[["human_symbol", "CHROM", "POS", "REF", "ALT", "impact", "cons_group", "hgvs_p", "AF"]].to_csv(
    f"{OUT}/mutations.tsv", sep="\t", index=False)

# --- parse GTF once: exons/CDS per transcript of interest ---
want_tx = {v[4]: k for k, v in GENES.items()}
exons = {t: [] for t in want_tx}
cds = {t: [] for t in want_tx}
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] not in ("exon", "CDS"):
            continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m:
            continue
        tx = m.group(1)
        if tx in want_tx:
            rec = (int(p[3]), int(p[4]))
            (exons if p[2] == "exon" else cds)[tx].append(rec)

# --- per gene: coverage + gene model ---
for gene, (chrom, gs, ge, strand, tx) in GENES.items():
    # coverage via samtools depth, binned
    region = f"{chrom}:{gs}-{ge}"
    try:
        out = subprocess.run([SAM, "depth", "-r", region, "-a", BAM],
                             capture_output=True, text=True, timeout=300).stdout
    except Exception as e:
        print(f"{gene}: depth FAILED {e}")
        out = ""
    pts = []
    for ln in out.splitlines():
        c, pos, d = ln.split("\t")
        pts.append((int(pos), int(d)))
    cov = pd.DataFrame(pts, columns=["pos", "depth"])
    # bin to ~150 bp
    if len(cov):
        span = ge - gs
        binw = max(50, int(span / 700))
        cov["bin"] = (cov["pos"] // binw) * binw + binw // 2
        covb = cov.groupby("bin")["depth"].mean().reset_index()
        covb.columns = ["pos", "depth"]
    else:
        covb = pd.DataFrame(columns=["pos", "depth"])
        binw = 0
    covb.to_csv(f"{OUT}/cov_{gene}.tsv", sep="\t", index=False)

    # gene model
    ex = pd.DataFrame(sorted(exons[tx]), columns=["start", "end"]) if exons[tx] else pd.DataFrame(columns=["start", "end"])
    cd = pd.DataFrame(sorted(cds[tx]), columns=["start", "end"]) if cds[tx] else pd.DataFrame(columns=["start", "end"])
    ex.to_csv(f"{OUT}/exons_{gene}.tsv", sep="\t", index=False)
    cd.to_csv(f"{OUT}/cds_{gene}.tsv", sep="\t", index=False)
    print(f"{gene}: cov_bins={len(covb)} (binw={binw}) exons={len(ex)} cds={len(cd)} strand={strand} tx={tx}")

print("done")
