#!/usr/bin/env Rscript
# Analysis B2: integrate promoter 5mC with copy-number state.
# Mechanism call per cancer gene:
#   homozygous deletion (log2r very low) | het-loss + promoter hyper (double hit) |
#   promoter hyper only (epigenetic silencing, CNV neutral) | amp + promoter hypo (activation)
suppressMessages({ library(data.table) })
WD <- "/workspace/hede_followup"
OUTD <- file.path(WD, "methyl"); dir.create(OUTD, showWarnings = FALSE)

prom   <- fread("/mnt/shared-workspace/shared/methyl/HeDe_promoter_5mC.tsv")  # gene,chrom(NC),pstart,pend,strand,mean5mC,n_cpg
cnv    <- fread("/mnt/shared-workspace/shared/HeDe_cnv_segments.tsv")          # ID,chr,start,end,nmark,log2r,width,state
cancer <- fread(file.path(WD, "ref/cancer_genes_rat.tsv"))

nc2chr <- function(nc) {
  num <- suppressWarnings(as.integer(sub("^NC_0*", "", nc)))
  out <- ifelse(!is.na(num) & num >= 86019 & num <= 86038, paste0("chr", num - 86018),
         ifelse(num == 86039, "chrX", ifelse(num == 86040, "chrY", NA_character_)))
  out
}
prom[, chr := nc2chr(chrom)]
prom <- prom[!is.na(chr)]
prom[, mid := as.integer((pstart + pend) / 2)]

# genome-wide promoter methylation baseline (for aberrant calling)
prom_base <- median(prom$mean5mC, na.rm = TRUE)
prom_hi   <- quantile(prom$mean5mC, 0.95, na.rm = TRUE)   # aberrant hyper threshold
prom_lo   <- quantile(prom$mean5mC, 0.05, na.rm = TRUE)   # aberrant hypo threshold
cat(sprintf("promoter 5mC: median=%.3f  q05=%.3f  q95=%.3f  (n=%d promoters)\n",
            prom_base, prom_lo, prom_hi, nrow(prom)))

# overlap promoter midpoint with CNV segments
cnv[, `:=`(start = as.integer(start), end = as.integer(end))]
setkey(cnv, chr, start, end)
prom_i <- prom[, .(chr, start = mid, end = mid, gene, mean5mC, n_cpg, pstart, pend)]
setkey(prom_i, chr, start, end)
ov <- foverlaps(prom_i, cnv[, .(chr, start, end, log2r, state)],
                by.x = c("chr","start","end"), by.y = c("chr","start","end"), nomatch = NULL)
ov <- unique(ov, by = c("gene","pstart","pend"))   # one segment per promoter (midpoint)

# aberrant-methylation flags (single-sample framework)
ov[, meth_flag := fcase(
  mean5mC >= prom_hi, "hyper",
  mean5mC <= prom_lo, "hypo",
  default = "typical")]
ov[, cn_state := state]

# mechanism classification
ov[, mechanism := fcase(
  log2r <= -2,                          "homozygous_loss",
  log2r <= -0.4 & meth_flag == "hyper", "het_loss+hyper (double-hit)",
  log2r <= -0.4,                        "het_loss",
  log2r >=  0.4 & meth_flag == "hypo",  "amp+hypo (activation)",
  log2r >=  0.4,                        "amplification",
  meth_flag == "hyper",                 "promoter_hyper_only",
  meth_flag == "hypo",                  "promoter_hypo_only",
  default = "none")]

fwrite(ov, file.path(OUTD, "promoter_5mC_x_cnv_all.tsv"), sep = "\t")
cat("\nmechanism distribution (all promoters):\n"); print(ov[, .N, by = mechanism][order(-N)])

# cancer-gene view
cg <- merge(ov, cancer, by.x = "gene", by.y = "rat_symbol")
setorder(cg, -is_hcc, role, human_symbol)
fwrite(cg, file.path(OUTD, "promoter_5mC_x_cnv_cancer.tsv"), sep = "\t")
cat("\ncancer genes with an aberrant mechanism (not 'none'):\n")
ab <- cg[mechanism != "none"]
print(ab[, .(human_symbol, gene, role, is_hcc, mean5mC = round(mean5mC,3),
             meth_flag, log2r, cn_state, mechanism)], nrows = 60)
cat("\ncancer genes by mechanism:\n"); print(cg[, .N, by = mechanism][order(-N)])

# key TSGs of interest
cat("\nkey loci:\n")
print(ov[gene %in% c("Cdkn2a","Cdkn2b","Tp53","Rb1","Pten","Axin1","Arid1a","Ctnnb1","Nfe2l2","Keap1"),
         .(gene, mean5mC = round(mean5mC,3), meth_flag, log2r, cn_state, mechanism)])
cat("DONE\n")
