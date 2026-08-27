#!/usr/bin/env bash
# Extract novel (ID==".") PASS SNV+indels from the dbSNP-annotated Clair3 VCF,
# then annotate functional consequences with the custom GRCr8 SnpEff database.
set -euo pipefail
export PATH=/workspace/.conda/envs/persist/bin:$PATH
JAVA=/workspace/.micromamba/envs/snpeff/bin/java
JAR=/workspace/.micromamba/envs/snpeff/share/snpeff-5.4.0c-0/snpEff.jar
CFG=/workspace/hede_followup/snpEff_grcr8.config
DATA=/workspace/hede_followup/snpeff_data
SRC=/mnt/shared-workspace/shared/variants/HeDe_annotated.vcf.gz
WD=/workspace/hede_followup/variants
mkdir -p "$WD"

echo "[1] Extract novel PASS variants (ID=='.' & FILTER=='PASS')"
zcat "$SRC" | awk '/^#/ {print; next} $3=="." && $7=="PASS"' > "$WD/novel_PASS.vcf"
echo "    novel PASS records: $(grep -vc '^#' "$WD/novel_PASS.vcf")"

echo "[2] Split SNV vs indel counts"
awk '!/^#/ { if (length($4)==1 && length($5)==1) s++; else i++ } END{print "    SNV="s, "indel="i}' "$WD/novel_PASS.vcf"

echo "[3] bgzip + index"
bgzip -f "$WD/novel_PASS.vcf"
tabix -f -p vcf "$WD/novel_PASS.vcf.gz"

echo "[4] SnpEff annotate (GRCr8 custom DB)"
"$JAVA" -Xmx8g -jar "$JAR" ann -v -c "$CFG" -dataDir "$DATA" GRCr8 "$WD/novel_PASS.vcf.gz" \
  2> "$WD/snpeff_ann.log" | bgzip > "$WD/novel_PASS.ann.vcf.gz"
echo "    annotated records: $(zcat "$WD/novel_PASS.ann.vcf.gz" | grep -vc '^#')"
echo "[5] ANN field present?"
zcat "$WD/novel_PASS.ann.vcf.gz" | grep -v '^#' | head -1 | cut -f8 | grep -o "ANN" || echo "NO ANN FIELD"
echo "DONE"
