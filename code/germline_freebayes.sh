#!/usr/bin/env bash
# Phase 3a: germline variant calling on the healthy-liver WGS with FreeBayes,
# parallelised per chromosome (single-threaded freebayes per chromosome, NJOBS
# concurrent), then concatenated and quality-filtered.
# Usage: germline_freebayes.sh <liver.markdup.bam> <outdir> [NJOBS] [CHROM_FILE]
# CHROM_FILE (optional): one accession per line; runs only those, skips final
# concat/filter (done centrally). Per-chr VCFs checkpoint to shared.
set -uo pipefail
BAM="$1"; OUT="$2"; NJOBS="${3:-8}"; CHROM_FILE="${4:-}"
REF=/workspace/ref/rn_GRCr8.fa
for E in wgs call s; do [ -d /workspace/.micromamba/envs/$E/bin ] && export PATH=/workspace/.micromamba/envs/$E/bin:$PATH; done
export PATH=/opt/conda/bin:$PATH
CKPT=/mnt/shared-workspace/shared/germline/per_chr
mkdir -p "$OUT/per_chr" /workspace/wgs/logs
LOG=/workspace/wgs/logs/germline_freebayes.log
echo "=== germline_freebayes $(date -u +%FT%TZ) ===" | tee "$LOG"

# chromosome accessions from the fasta index (NC_ only, i.e. placed chromosomes)
if [ -n "$CHROM_FILE" ]; then CHROMS=$(cat "$CHROM_FILE"); else
CHROMS=$(awk '$1 ~ /^NC_/ {print $1}' "$REF.fai"); fi
echo "chromosomes: $(echo "$CHROMS" | wc -l)" | tee -a "$LOG"

run_chr() {
  c="$1"
  mkdir -p "$CKPT"
  if [ -s "$CKPT/$c.vcf.gz" ] && [ -s "$CKPT/$c.vcf.gz.tbi" ]; then
    echo "  $c already checkpointed, skip" >>"$LOG"; return 0
  fi
  freebayes -f "$REF" --region "$c" "$BAM" 2>>"$LOG" | \
    bcftools view -Oz -o "$OUT/per_chr/$c.vcf.gz" 2>>"$LOG"
  bcftools index -t "$OUT/per_chr/$c.vcf.gz" 2>>"$LOG"
  cp "$OUT/per_chr/$c.vcf.gz" "$OUT/per_chr/$c.vcf.gz.tbi" "$CKPT/"
  echo "  $c done: $(bcftools view -H "$OUT/per_chr/$c.vcf.gz" | wc -l) records" >>"$LOG"
}
export -f run_chr; export REF BAM OUT LOG CKPT

echo "$CHROMS" | xargs -P "$NJOBS" -I{} bash -c 'run_chr "$@"' _ {}

if [ -n "$CHROM_FILE" ]; then echo "subset mode: skipping concat" | tee -a "$LOG"; exit 0; fi
bcftools concat -Oz -o "$OUT/HeDe_liver_germline_raw.vcf.gz" \
  $(for c in $CHROMS; do echo "$OUT/per_chr/$c.vcf.gz"; done) >>"$LOG" 2>&1
bcftools index -t "$OUT/HeDe_liver_germline_raw.vcf.gz" >>"$LOG" 2>&1

# standard hard filter (documented; no truth set exists for F344 rat)
bcftools filter -Oz -s LOWQUAL -e 'QUAL < 30' \
  "$OUT/HeDe_liver_germline_raw.vcf.gz" \
  -o "$OUT/HeDe_liver_germline.vcf.gz" >>"$LOG" 2>&1
bcftools index -t "$OUT/HeDe_liver_germline.vcf.gz" >>"$LOG" 2>&1

echo "raw records:  $(bcftools view -H "$OUT/HeDe_liver_germline_raw.vcf.gz" | wc -l)" | tee -a "$LOG"
echo "PASS records: $(bcftools view -H -f PASS "$OUT/HeDe_liver_germline.vcf.gz" | wc -l)" | tee -a "$LOG"
echo "=== DONE $(date -u +%FT%TZ) ===" | tee -a "$LOG"
