#!/bin/bash
# somatic_mutect2.sh — Phase 4 somatic SNV/indel via Mutect2 (GATK 4.6), tumor H_1 vs
# strain-matched liver RL_1. Replaces Strelka2 (2.9.10 segfaults on minimap2 BAMs).
# Scatter Mutect2 per chromosome -> GatherVcfs -> FilterMutectCalls -> PASS VCFs.
# Usage: somatic_mutect2.sh <tumor.bam> <normal.bam> <outdir> <threads>
set -uo pipefail
TUM="$1"; NOR="$2"; OUT="$3"; THREADS="${4:-6}"
REF=/workspace/ref/rn_GRCr8.fa
export PATH=/workspace/.micromamba/envs/gatk46/bin:/workspace/.micromamba/envs/som/bin:/opt/conda/bin:$PATH
export LD_LIBRARY_PATH=/workspace/.micromamba/envs/gatk46/lib:${LD_LIBRARY_PATH:-}
mkdir -p "$OUT/mutect2" /workspace/wgs/logs
LOG=/workspace/wgs/logs/somatic_mutect2.log
echo "=== somatic_mutect2 $(date -u +%FT%TZ) ===" | tee "$LOG"

# chromosome list (22 NC_ accessions, .fai order)
mapfile -t CHROMS < <(grep '^NC_' "${REF}.fai" | cut -f1)
echo "chromosomes: ${#CHROMS[@]}" | tee -a "$LOG"

# scatter Mutect2 per chromosome (each single-threaded region; run THREADS in parallel)
run_chr(){
  local C="$1"
  gatk Mutect2 -R "$REF" \
    -I "$TUM" -I "$NOR" -normal RL_1 \
    -L "$C" \
    --native-pair-hmm-threads 1 \
    -O "$OUT/mutect2/$C.vcf.gz" >>"$LOG" 2>&1 \
    && echo "  done $C" >>"$LOG" || echo "  FAIL $C" >>"$LOG"
}
export -f run_chr; export REF TUM NOR OUT LOG
printf '%s\n' "${CHROMS[@]}" | xargs -P "$THREADS" -I{} bash -c 'run_chr "$@"' _ {}
echo "per-chr VCFs: $(ls "$OUT"/mutect2/NC_*.vcf.gz 2>/dev/null | wc -l) / ${#CHROMS[@]}" | tee -a "$LOG"

# gather
gatk GatherVcfs $(for C in "${CHROMS[@]}"; do echo -n "-I $OUT/mutect2/$C.vcf.gz "; done) \
  -O "$OUT/mutect2/HeDe_WGS_mutect2_unfiltered.vcf.gz" >>"$LOG" 2>&1 || { echo "ABORT gather" | tee -a "$LOG"; exit 1; }

# merge per-chr stats + f1r2 (none here) then filter
gatk MergeMutectStats $(for C in "${CHROMS[@]}"; do echo -n "--stats $OUT/mutect2/$C.vcf.gz.stats "; done) \
  -O "$OUT/mutect2/merged.stats" >>"$LOG" 2>&1
gatk FilterMutectCalls -R "$REF" \
  -V "$OUT/mutect2/HeDe_WGS_mutect2_unfiltered.vcf.gz" \
  --stats "$OUT/mutect2/merged.stats" \
  -O "$OUT/HeDe_WGS_somatic_mutect2.vcf.gz" >>"$LOG" 2>&1 || { echo "ABORT filter" | tee -a "$LOG"; exit 1; }

V="$OUT/HeDe_WGS_somatic_mutect2.vcf.gz"
echo "mutect2 records (all):  $(bcftools view -H "$V" | wc -l)" | tee -a "$LOG"
echo "mutect2 records (PASS): $(bcftools view -H -f PASS "$V" | wc -l)" | tee -a "$LOG"
echo "  SNVs PASS:   $(bcftools view -H -f PASS -v snps "$V" | wc -l)" | tee -a "$LOG"
echo "  indels PASS: $(bcftools view -H -f PASS -v indels "$V" | wc -l)" | tee -a "$LOG"

# positive control: Tp53 p.Arg271Ser = NC_086028.1:54807961 C>A
echo "--- Tp53 positive control (NC_086028.1:54807961):" | tee -a "$LOG"
bcftools view -H -r NC_086028.1:54807961 "$V" 2>/dev/null | tee -a "$LOG" || echo "  NOT CALLED" | tee -a "$LOG"
echo "=== DONE $(date -u +%FT%TZ) ===" | tee -a "$LOG"
