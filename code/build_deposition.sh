#!/usr/bin/env bash
# Build the deposition package for the HeDe PacBio Data Descriptor.
# Organizes depositable files into zenodo/ and sra/ bundles with README, manifest, checksums.
set -euo pipefail

SW=/mnt/shared-workspace/shared
FU=/workspace/hede_followup
RES=/mnt/results
BUILD=/workspace/deposition
Z=$BUILD/zenodo/HeDe_PacBio_GRCr8_dataset
E=$BUILD/sra

rm -rf "$BUILD"
mkdir -p "$Z"/{assembly,variants,sv,cnv,methylation,tables,code} "$E"

echo "[1/6] assembly (gzip)"
gzip -c "$SW/compleasm_input/HeDe_assembly.fasta" > "$Z/assembly/HeDe_assembly.fasta.gz"

echo "[2/6] variants + SV + CNV"
cp "$FU/variants/novel_PASS.vcf.gz"      "$Z/variants/HeDe_novel_PASS.vcf.gz"
cp "$FU/variants/novel_PASS.vcf.gz.tbi"  "$Z/variants/HeDe_novel_PASS.vcf.gz.tbi"
cp "$FU/variants/novel_PASS.ann.vcf.gz"  "$Z/variants/HeDe_novel_PASS.ann.vcf.gz"
cp "$SW/sv/HeDe_sniffles2.vcf.gz"        "$Z/sv/"
cp "$SW/sv/HeDe_sniffles2.vcf.gz.tbi"    "$Z/sv/"
cp "$SW/sv/HeDe_sv_table.tsv"            "$Z/sv/"
cp "$SW/HeDe_cnv_segments.tsv"           "$Z/cnv/"
cp "$SW/HeDe_cnv_20kb_bins.tsv"          "$Z/cnv/"

echo "[3/6] methylation"
cp "$SW/methyl/HeDe_cpg_5mC.bed.gz"      "$Z/methylation/"
cp "$SW/methyl/HeDe_promoter_5mC.tsv"    "$Z/methylation/"
cp "$RES/followup_tables/cgi_5mC.tsv"           "$Z/methylation/"
cp "$RES/followup_tables/retrotransposon_5mC.tsv" "$Z/methylation/"
cp "$RES/followup_tables/imprinted_5mC.tsv"     "$Z/methylation/"

echo "[4/6] processed tables"
cp "$RES/followup_tables/cancer_mutations_functional_impact.tsv" "$Z/tables/"
cp "$RES/followup_tables/cancer_coding_mutations.tsv"            "$Z/tables/"
cp "$RES/followup_tables/coding_mutations_all.tsv"               "$Z/tables/"
cp "$RES/followup_tables/promoter_5mC_x_cnv_cancer.tsv"          "$Z/tables/"
cp "$RES/mutational_spectrum_6class.csv"  "$Z/tables/"
cp "$RES/mutational_spectrum_96.csv"      "$Z/tables/"
cp "$RES/signature_contributions_targeted.csv" "$Z/tables/"
cp "$SW/variants/HeDe_signature_contributions.csv" "$Z/tables/"

echo "[5/6] code"
cp "$FU"/*.R "$FU"/*.py "$FU"/*.sh "$Z/code/" 2>/dev/null || true

echo "[6/6] manifest + checksums"
cd "$Z"
find . -type f ! -name "MD5SUMS" ! -name "MANIFEST.tsv" ! -name "README.md" -print0 | sort -z | xargs -0 md5sum > MD5SUMS
# manifest: path, size_bytes, md5
python3 - <<'PY'
import os, hashlib
root="."
rows=[]
for dp,_,fns in os.walk(root):
    for fn in fns:
        if fn in ("MD5SUMS","MANIFEST.tsv","README.md"): continue
        p=os.path.join(dp,fn)
        rel=os.path.relpath(p,root)
        rows.append((rel, os.path.getsize(p)))
rows.sort()
with open("MANIFEST.tsv","w") as f:
    f.write("file\tsize_bytes\n")
    for rel,sz in rows:
        f.write(f"{rel}\t{sz}\n")
print(f"manifest rows: {len(rows)}")
PY
echo "BUILD OK"
du -sh "$Z"
