#!/usr/bin/env bash
# Phase 4: somatic SNV/indel calling, tumor(H_1) vs strain-matched liver(RL_1).
# Manta (SV candidates for Strelka2) -> Strelka2 somatic -> PASS SNV/indel VCFs.
# Usage: somatic_manta_strelka.sh <tumor.bam> <normal.bam> <outdir> <threads>
set -uo pipefail
TUM="$1"; NOR="$2"; OUT="$3"; THREADS="${4:-14}"
REF=/workspace/ref/rn_GRCr8.fa
ENVROOT=/workspace/.micromamba/envs
MENV=""
for E in som wgs; do [ -x "$ENVROOT/$E/bin/configManta.py" ] && MENV="$ENVROOT/$E" && break; done
[ -z "$MENV" ] && { echo "ABORT: no env with configManta.py" | tee -a /dev/stderr; exit 1; }
export PATH=$MENV/bin:/opt/conda/bin:$PATH
export LD_LIBRARY_PATH=$MENV/lib:${LD_LIBRARY_PATH:-}
mkdir -p "$OUT/manta" "$OUT/strelka" /workspace/wgs/logs
LOG=/workspace/wgs/logs/somatic.log
echo "=== somatic_manta_strelka $(date -u +%FT%TZ) ===" | tee "$LOG"

# Manta: somatic SV (also produces candidate indels for Strelka2)
configManta.py --tumorBam "$TUM" --normalBam "$NOR" --referenceFasta "$REF" \
  --runDir "$OUT/manta" >>"$LOG" 2>&1 || { echo "ABORT: manta config" | tee -a "$LOG"; exit 1; }
"$OUT/manta/runWorkflow.py" -m local -j "$THREADS" >>"$LOG" 2>&1 || { echo "ABORT: manta run" | tee -a "$LOG"; exit 1; }
echo "manta done: $(bcftools view -H "$OUT/manta/results/variants/somaticSV.vcf.gz" | wc -l) somatic SV records" | tee -a "$LOG"

# Strelka2: somatic SNV+indel, guided by Manta candidate indels
configureStrelkaSomaticWorkflow.py --tumorBam "$TUM" --normalBam "$NOR" \
  --referenceFasta "$REF" \
  --indelCandidates "$OUT/manta/results/variants/candidateSmallIndels.vcf.gz" \
  --runDir "$OUT/strelka" >>"$LOG" 2>&1 || { echo "ABORT: strelka config" | tee -a "$LOG"; exit 1; }
"$OUT/strelka/runWorkflow.py" -m local -j "$THREADS" >>"$LOG" 2>&1 || { echo "ABORT: strelka run" | tee -a "$LOG"; exit 1; }

SNV="$OUT/strelka/results/variants/somatic.snvs.vcf.gz"
IND="$OUT/strelka/results/variants/somatic.indels.vcf.gz"
echo "strelka SNVs (all):    $(bcftools view -H "$SNV" | wc -l)" | tee -a "$LOG"
echo "strelka SNVs (PASS):   $(bcftools view -H -f PASS "$SNV" | wc -l)" | tee -a "$LOG"
echo "strelka indels (all):  $(bcftools view -H "$IND" | wc -l)" | tee -a "$LOG"
echo "strelka indels (PASS): $(bcftools view -H -f PASS "$IND" | wc -l)" | tee -a "$LOG"

# positive control: Tp53 p.Arg271Ser = NC_086028.1:54807961 C>A must be present
echo "--- Tp53 positive control (NC_086028.1:54807961 C>A):" | tee -a "$LOG"
bcftools view -H -r NC_086028.1:54807961 "$SNV" 2>/dev/null | tee -a "$LOG" || true
bcftools view -H -r NC_086028.1:54807960-54807962 "$IND" 2>/dev/null | tee -a "$LOG" || true

cp "$SNV" "$OUT/HeDe_WGS_somatic_snvs.vcf.gz"
cp "$IND" "$OUT/HeDe_WGS_somatic_indels.vcf.gz"
cp "$OUT/manta/results/variants/somaticSV.vcf.gz" "$OUT/HeDe_WGS_somatic_sv_manta.vcf.gz"
echo "=== DONE $(date -u +%FT%TZ) ===" | tee -a "$LOG"
