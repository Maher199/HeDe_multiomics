#!/bin/bash
# align_local_minimap2.sh <SAMPLE_LANE> <R1.fq.gz> <R2.fq.gz>
# LOCAL fallback alignment: minimap2 -ax sr at SAFE thread count (-t 6, ~37% load)
# piped straight to samtools sort (-@ 2 -m 1G) -> sorted BAM + index + flagstat.
# Checkpoint-safe: per-lane output; rerun only missing lanes.
set -euo pipefail
export MAMBA_ROOT_PREFIX=/workspace/.micromamba
eval "$(/workspace/.micromamba/bin/micromamba shell hook -s bash 2>/dev/null)" || true
micromamba activate wgs   # provides samtools; minimap2 added below

TAG="$1"; R1="$2"; R2="$3"
REF=/workspace/ref/rn_GRCr8.fa
OUTDIR=/workspace/wgs/bam
mkdir -p "$OUTDIR"
OUT="$OUTDIR/${TAG}.sorted.bam"
MM2=/workspace/tools/minimap2-2.28/minimap2
RG="@RG\\tID:${TAG}\\tSM:${TAG%_L*}\\tLB:${TAG%_L*}\\tPL:ILLUMINA\\tPU:${TAG}"

if [ -s "$OUT" ] && [ -s "$OUT.bai" ]; then
  echo "[$(date -Is)] $OUT exists, skipping"; exit 0
fi
echo "[$(date -Is)] aligning $TAG (minimap2 -t 6 local)"
"$MM2" -ax sr -t 6 -R "$RG" "$REF" "$R1" "$R2" 2> "$OUTDIR/${TAG}.minimap2.log" \
  | samtools sort -@ 2 -m 1G -o "$OUT" -
samtools index "$OUT"
samtools flagstat "$OUT" > "$OUTDIR/${TAG}.flagstat.txt"
samtools quickcheck -v "$OUT" && echo "[$(date -Is)] DONE $TAG"
