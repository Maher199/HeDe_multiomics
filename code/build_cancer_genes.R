#!/usr/bin/env Rscript
# Build a rat (GRCr8) cancer-gene set from the intOGen Compendium of Cancer Genes.
suppressMessages({library(data.table)})
ref_dir  <- "/workspace/hede_followup/ref"
intogen  <- fread(file.path(ref_dir, "2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv"), fill = TRUE)
ratgenes <- fread(file.path(ref_dir, "GRCr8_genes.bed"),
                  col.names = c("chrom","start","end","strand","rat_symbol"))

# Per-gene summary across cohorts (human symbol)
setnames(intogen, c("SYMBOL","CANCER_TYPE","ROLE"), c("human_symbol","ctype","role"))
summ <- intogen[, .(
    n_cohorts = uniqueN(COHORT),
    n_ctypes  = uniqueN(ctype),
    ctypes    = paste(sort(unique(ctype)), collapse=","),
    is_hcc    = any(ctype == "HC"),
    role      = { r <- role[role %in% c("Act","LoF")]
                  if (length(unique(r))==1) unique(r) else "ambiguous" }
  ), by = human_symbol]

ratgenes[, usym := toupper(rat_symbol)]
summ[, usym := toupper(human_symbol)]
merged <- merge(summ, ratgenes, by = "usym", allow.cartesian = TRUE)
merged <- merged[!is.na(chrom)]
setcolorder(merged, c("human_symbol","rat_symbol","usym","chrom","start","end","strand","role","is_hcc","n_ctypes","ctypes"))
setorder(merged, chrom, start)

fwrite(merged, file.path(ref_dir, "cancer_genes_rat.tsv"), sep = "\t")
cat("intOGen unique driver genes:", nrow(summ), "\n")
cat("Mapped to rat genes (unique symbols):", uniqueN(merged$usym), "\n")
cat("HCC (HC) drivers mapped:", uniqueN(merged$usym[merged$is_hcc]), "\n")
cat("Role distribution (mapped gene rows):\n"); print(merged[, .N, by=role])
print(merged[usym %in% c("TP53","KRAS","HRAS","CDKN2A","MYC","CTNNB1","PTEN","RB1"),
             .(human_symbol,rat_symbol,chrom,start,end,role,is_hcc)])
