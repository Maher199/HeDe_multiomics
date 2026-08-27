#!/bin/bash
# germline_concat.sh — concat 22 per-chr FreeBayes VCFs in .fai order, QUAL>=30 filter.
set -euo pipefail
SW=/mnt/shared-workspace/shared
export MAMBA_ROOT_PREFIX=/workspace/.micromamba
log(){ echo "[$(date -Is)] [concat] $*"; }

# need bcftools: reuse 'call' env if present, else create
BCF=""
for E in call wgs som; do [ -x /workspace/.micromamba/envs/$E/bin/bcftools ] && BCF=/workspace/.micromamba/envs/$E/bin/bcftools && break; done
if [ -z "$BCF" ]; then
  log "creating call env (bcftools)"
  /usr/bin/micromamba create -y -n call -c conda-forge -c bioconda bcftools > /workspace/mm_call.log 2>&1
  BCF=/workspace/.micromamba/envs/call/bin/bcftools
fi
log "bcftools: $BCF"

mkdir -p /workspace/germline
# order = .fai NC_ accessions
mapfile -t CHROMS < <(grep '^NC_' "$SW/rn_GRCr8.fa.fai" | cut -f1)
log "concat ${#CHROMS[@]} chroms in .fai order"
FILES=(); for C in "${CHROMS[@]}"; do FILES+=("$SW/germline/per_chr/$C.vcf.gz"); done
"$BCF" concat -a -O z -o /workspace/germline/HeDe_liver_germline_raw.vcf.gz "${FILES[@]}"
"$BCF" index -t /workspace/germline/HeDe_liver_germline_raw.vcf.gz
RAW=$("$BCF" view -H /workspace/germline/HeDe_liver_germline_raw.vcf.gz | wc -l)
log "raw records: $RAW"

"$BCF" filter -s LOWQUAL -e 'QUAL<30' -O z -o /workspace/germline/HeDe_liver_germline.vcf.gz /workspace/germline/HeDe_liver_germline_raw.vcf.gz
"$BCF" index -t /workspace/germline/HeDe_liver_germline.vcf.gz
PASS=$("$BCF" view -H -f PASS /workspace/germline/HeDe_liver_germline.vcf.gz | wc -l)
log "PASS (QUAL>=30) records: $PASS"

# checkpoint
cp /workspace/germline/HeDe_liver_germline_raw.vcf.gz* "$SW/germline/"
cp /workspace/germline/HeDe_liver_germline.vcf.gz* "$SW/germline/"
log "checkpointed. DONE raw=$RAW pass=$PASS"
