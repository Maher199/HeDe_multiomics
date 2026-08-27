#!/usr/bin/env python3
"""Score all unique cancer-gene SNVs with rat SIFT via Ensembl VEP REST (GRCr8)."""
import json, time, urllib.request
import pandas as pd

WD = "/workspace/hede_followup"
cc = pd.read_csv(f"{WD}/variants/cancer_coding_mutations.tsv", sep="\t")
u = cc.drop_duplicates(subset=["human_symbol", "CHROM", "POS"]).copy()

# NC_* -> chrN  (NC_086019=chr1 ... NC_086038=chr20, NC_086039=chrX, NC_086040=chrY)
def nc2chr(nc):
    n = int(nc.split(".")[0].replace("NC_0860", ""))
    if n <= 38:
        return str(n - 18)
    return {39: "X", 40: "Y"}[n]

u["chr"] = u["CHROM"].map(nc2chr)
# SNVs only (single-base REF and ALT) get SIFT
snv = u[(u["REF"].str.len() == 1) & (u["ALT"].str.len() == 1)].copy()
print(f"unique cancer variants: {len(u)}; SNVs to score: {len(snv)}")

# Build VEP region-format strings: "chr start end ref/alt strand"
vep_ids = []
for _, r in snv.iterrows():
    vep_ids.append(f"{r['chr']} {r['POS']} {r['POS']} {r['REF']}/{r['ALT']} 1")

URL = "https://rest.ensembl.org/vep/rattus_norvegicus/region"
HDRS = {"Content-Type": "application/json", "Accept": "application/json"}

def post_batch(ids):
    body = json.dumps({"variants": ids}).encode()
    req = urllib.request.Request(URL, data=body, headers=HDRS)
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.load(resp)

# batch in chunks of 50
results = {}
B = 50
for i in range(0, len(vep_ids), B):
    chunk = vep_ids[i:i+B]
    for attempt in range(4):
        try:
            out = post_batch(chunk)
            break
        except Exception as e:
            print(f"  batch {i//B} attempt {attempt} error: {e}; retrying")
            time.sleep(2 * (attempt + 1))
    else:
        out = []
    for v in out:
        # pick the transcript consequence matching the gene with a missense/SIFT call
        best = None
        for tc in v.get("transcript_consequences", []):
            if tc.get("sift_score") is not None:
                if best is None or tc.get("canonical") == 1:
                    best = tc
        if best is not None:
            key = v.get("input", "")
            results[key] = (best.get("sift_prediction"), best.get("sift_score"),
                            best.get("gene_symbol"), best.get("amino_acids"))
    print(f"  scored {min(i+B, len(vep_ids))}/{len(vep_ids)}")
    time.sleep(0.4)

# map back
def lookup(r):
    key = f"{r['chr']} {r['POS']} {r['POS']} {r['REF']}/{r['ALT']} 1"
    return results.get(key, (None, None, None, None))

snv[["rat_sift_pred", "rat_sift_score", "sift_gene", "sift_aa"]] = snv.apply(
    lambda r: pd.Series(lookup(r)), axis=1)

out = snv[["human_symbol", "CHROM", "POS", "REF", "ALT", "rat_sift_pred", "rat_sift_score", "sift_gene", "sift_aa"]]
out.to_csv(f"{WD}/variants/rat_sift_cancer.tsv", sep="\t", index=False)
print("\nSIFT prediction counts:")
print(snv["rat_sift_pred"].value_counts(dropna=False))
print("\nwrote variants/rat_sift_cancer.tsv")
