#!/bin/bash
# merge_markdup2.sh <SAMPLE> <TAG1> [TAG2 ...]
# Merge checkpointed unit BAMs from shared/bam/, markdup, index, stats, mosdepth.
# Uses the lightweight 's' env (samtools) + installs mosdepth into it.
set -euo pipefail
SAMPLE="$1"; shift
SW=/mnt/shared-workspace/shared
log(){ echo "[$(date -Is)] [$SAMPLE] $*"; }

export MAMBA_ROOT_PREFIX=/workspace/.micromamba
if [ ! -x /workspace/.micromamba/envs/s/bin/samtools ]; then
  /usr/bin/micromamba create -y -n s -c conda-forge -c bioconda samtools > /workspace/mm_create.log 2>&1
fi
if [ ! -x /workspace/.micromamba/envs/s/bin/mosdepth ]; then
  /usr/bin/micromamba install -y -n s -c conda-forge -c bioconda mosdepth >> /workspace/mm_create.log 2>&1
fi
eval "$(/usr/bin/micromamba shell hook -s bash)"; micromamba activate s

mkdir -p /workspace/post
cd /workspace/post
# localize inputs
IN=()
for T in "$@"; do
  [ -s "${T}.sorted.bam" ] || { log "fetch $T"; cp "$SW/bam/${T}.sorted.bam" .; cp "$SW/bam/${T}.sorted.bam.bai" .; }
  IN+=("${T}.sorted.bam")
done

log "merge ${#IN[@]} units"
rm -f ${SAMPLE}.merged.bam
if [ "${#IN[@]}" -eq 1 ]; then cp "${IN[0]}" ${SAMPLE}.merged.bam; else samtools merge -f -@ 4 -o ${SAMPLE}.merged.bam "${IN[@]}"; fi
log "namesort for fixmate (minimap2 lacks MC tag; markdup needs it)"
samtools sort -n -@ 4 -m 1G -o ${SAMPLE}.namesort.bam ${SAMPLE}.merged.bam
rm -f ${SAMPLE}.merged.bam
log "fixmate"
samtools fixmate -@ 4 -m ${SAMPLE}.namesort.bam ${SAMPLE}.fixmate.bam
rm -f ${SAMPLE}.namesort.bam
log "coordinate re-sort"
samtools sort -@ 4 -m 1G -o ${SAMPLE}.fixmate.sorted.bam ${SAMPLE}.fixmate.bam
rm -f ${SAMPLE}.fixmate.bam
log "markdup"
samtools markdup -@ 4 ${SAMPLE}.fixmate.sorted.bam ${SAMPLE}.markdup.bam
samtools index ${SAMPLE}.markdup.bam
rm -f ${SAMPLE}.fixmate.sorted.bam
log "stats"
samtools flagstat ${SAMPLE}.markdup.bam > ${SAMPLE}.markdup.flagstat.txt
samtools stats ${SAMPLE}.markdup.bam > ${SAMPLE}.markdup.stats.txt
log "mosdepth"
mosdepth -t 4 -n -b 20000 ${SAMPLE}.20kb ${SAMPLE}.markdup.bam
samtools quickcheck -v ${SAMPLE}.markdup.bam
log "checkpoint"
mkdir -p "$SW/bam/final"
cp ${SAMPLE}.markdup.bam ${SAMPLE}.markdup.bam.bai ${SAMPLE}.markdup.flagstat.txt ${SAMPLE}.markdup.stats.txt ${SAMPLE}.20kb.regions.bed.gz ${SAMPLE}.20kb.mosdepth.*.txt "$SW/bam/final/" 2>/dev/null || true
ls -la "$SW/bam/final/"
log "DONE"
