#!/usr/bin/env python3
"""Extract per-gene data for fig8 browser views (second-tier HCC drivers): coverage, gene model, mutations."""
import subprocess, gzip, re, os
import pandas as pd

WD = "/workspace/hede_followup"
BAM = "/mnt/shared-workspace/shared/HeDe_on_GRCr8.sorted.bam"
GTF = "/workspace/hede_followup/ref/GCF_036323735.1_GRCr8_genomic.gtf.gz"
SAM = "/workspace/.conda/envs/persist/bin/samtools"
OUT = os.path.join(WD, "fig8_data")
os.makedirs(OUT, exist_ok=True)

# gene -> (chrom, start, end, strand, transcript)  [ordered by intOGen n_ctypes]
GENES = {
    "KMT2C": ("NC_086022.1", 10353698, 10755965, "+", "XM_039108963.2"),
    "ARID1B":("NC_086019.1", 47973199, 48328793, "+", "NM_001419802.1"),
    "TET1":  ("NC_086038.1", 25766806, 25839598, "+", "XM_039099324.2"),
    "CLTC":  ("NC_086028.1", 72014984, 72073308, "-", "XM_063269683.1"),
    "WNK2":  ("NC_086035.1", 15709648, 15818874, "+", "XM_006253715.5"),
    "DHX9":  ("NC_086031.1", 68152813, 68189580, "-", "NM_001395556.1"),
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

    ex = pd.DataFrame(sorted(exons[tx]), columns=["start", "end"]) if exons[tx] else pd.DataFrame(columns=["start", "end"])
    cd = pd.DataFrame(sorted(cds[tx]), columns=["start", "end"]) if cds[tx] else pd.DataFrame(columns=["start", "end"])
    ex.to_csv(f"{OUT}/exons_{gene}.tsv", sep="\t", index=False)
    cd.to_csv(f"{OUT}/cds_{gene}.tsv", sep="\t", index=False)
    print(f"{gene}: cov_bins={len(covb)} (binw={binw}) exons={len(ex)} cds={len(cd)} strand={strand} tx={tx}")

print("done")
