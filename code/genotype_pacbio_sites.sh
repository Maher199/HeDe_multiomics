#!/usr/bin/env bash
# Phase 3b: targeted genotyping of the 2,234,197 PacBio novel-SNV sites in both
# WGS samples (HeDe tumor + healthy liver) via bcftools mpileup + call -m.
# Produces per-sample GT at every PacBio site for the orthogonal-validation table.
# Usage: genotype_pacbio_sites.sh <tumor.bam> <liver.bam> <outdir>
set -uo pipefail
TUM="$1"; NOR="$2"; OUT="$3"
REF=/workspace/ref/rn_GRCr8.fa
SITES=/workspace/wgs/pacbio_novel_sites.nohdr.tsv   # header-stripped CHROM POS REF ALT
for E in wgs call s; do [ -d /workspace/.micromamba/envs/$E/bin ] && export PATH=/workspace/.micromamba/envs/$E/bin:$PATH; done
export PATH=/opt/conda/bin:$PATH
mkdir -p "$OUT" /workspace/wgs/logs
LOG=/workspace/wgs/logs/genotype_pacbio_sites.log
echo "=== genotype_pacbio_sites $(date -u +%FT%TZ) ===" | tee "$LOG"

# strip header line if present (mpileup -T cannot parse it)
if [ ! -s "$SITES" ]; then
  SRC=/workspace/wgs/pacbio_novel_sites.tsv
  [ -s "$SRC" ] || { echo "ABORT: $SRC missing" | tee -a "$LOG"; exit 1; }
  awk 'NR==1 && $1=="chr"{next} {print}' "$SRC" > "$SITES"
fi

# mpileup at exactly the PacBio sites, both samples together; -a FORMAT/AD for allele depths
bcftools mpileup -f "$REF" -T "$SITES" -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 \
  "$TUM" "$NOR" 2>>"$LOG" | \
  bcftools call -m -Oz -o "$OUT/HeDe_WGS_at_pacbio_sites.vcf.gz" 2>>"$LOG"
bcftools index -t "$OUT/HeDe_WGS_at_pacbio_sites.vcf.gz" >>"$LOG" 2>&1

echo "sites genotyped: $(bcftools view -H "$OUT/HeDe_WGS_at_pacbio_sites.vcf.gz" | wc -l) / 2234197" | tee -a "$LOG"
echo "=== DONE $(date -u +%FT%TZ) ===" | tee -a "$LOG"
