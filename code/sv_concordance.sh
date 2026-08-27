#!/usr/bin/env bash
# Phase 6: SV concordance between WGS (Manta, tumor/normal) and PacBio (Sniffles2,
# tumor-only) via SURVIVOR merge: type-matched, strand-agnostic, max breakpoint
# distance 500 bp, min SV size 50 bp. SUPP_VEC then classifies each merged SV as
# Manta-only / Sniffles-only / both.
# Usage: sv_concordance.sh <manta_somaticSV.vcf.gz> <sniffles.vcf.gz> <outdir>
set -uo pipefail
MANTA="$1"; SNIF="$2"; OUT="$3"
export PATH=/workspace/.micromamba/envs/wgs/bin:/opt/conda/bin:$PATH
mkdir -p "$OUT" /workspace/wgs/logs
LOG=/workspace/wgs/logs/sv_concordance.log
echo "=== sv_concordance $(date -u +%FT%TZ) ===" | tee "$LOG"

# SURVIVOR wants uncompressed VCFs; restrict to PASS records
bcftools view -f PASS -Ov -o "$OUT/manta_pass.vcf" "$MANTA" 2>>"$LOG"
bcftools view -f PASS -Ov -o "$OUT/sniffles_pass.vcf" "$SNIF" 2>>"$LOG"
echo "manta PASS:    $(grep -vc '^#' "$OUT/manta_pass.vcf")" | tee -a "$LOG"
echo "sniffles PASS: $(grep -vc '^#' "$OUT/sniffles_pass.vcf")" | tee -a "$LOG"

printf '%s\n%s\n' "$OUT/manta_pass.vcf" "$OUT/sniffles_pass.vcf" > "$OUT/vcf_list.txt"
# merge: maxdist 500, min callers 1, type-matched 1, strand-agnostic 0, est-dist 0, minsize 50
SURVIVOR merge "$OUT/vcf_list.txt" 500 1 1 0 0 50 "$OUT/merged_survivor.vcf" >>"$LOG" 2>&1
echo "merged SVs: $(grep -vc '^#' "$OUT/merged_survivor.vcf")" | tee -a "$LOG"

/opt/conda/bin/Rscript /workspace/wgs/sv_concordance.R "$OUT/merged_survivor.vcf" "$OUT" | tee -a "$LOG"
echo "=== DONE $(date -u +%FT%TZ) ===" | tee -a "$LOG"
