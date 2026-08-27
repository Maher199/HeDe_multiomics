#!/usr/bin/env Rscript
# Consolidate rat SIFT + human-ortholog scores + CNV/LOH state into one table.
suppressMessages(library(data.table))
WD <- "/workspace/hede_followup"

cc <- fread(file.path(WD, "variants/cancer_coding_mutations.tsv"))
u <- unique(cc, by = c("human_symbol", "CHROM", "POS"))

# NC_* -> chrN
nc2chr <- function(nc) {
  n <- as.integer(sub("NC_0860", "", sub("\\..*", "", nc)))
  ifelse(n <= 38, as.character(n - 18), ifelse(n == 39, "X", "Y"))
}
u[, chr := paste0("chr", nc2chr(CHROM))]  # CNV file uses chrN naming

# --- rat SIFT ---
sift <- fread(file.path(WD, "variants/rat_sift_cancer.tsv"))
sift <- sift[, .(human_symbol, CHROM, POS, rat_sift_pred, rat_sift_score)]
u <- merge(u, sift, by = c("human_symbol", "CHROM", "POS"), all.x = TRUE)

# --- human ortholog scores (driver set + deep-dive) ---
ho <- fread(file.path(WD, "variants/human_ortholog_scores.tsv"))
dd <- fread(file.path(WD, "variants/deepdive_human_scores.tsv"))
dd <- dd[, .(gene, rat_hgvsp, human_hgvsp, cadd_phred, am_class, am_score, human_sift, status)]
ho <- ho[, .(gene, rat_hgvsp, human_hgvsp, cadd_phred, am_class, am_score, human_sift, status)]
hum <- rbind(ho, dd, fill = TRUE)
hum <- hum[!is.na(human_hgvsp) | status != "OK" | TRUE]  # keep all
# collapse to one row per gene+rat_hgvsp
hum <- unique(hum, by = c("gene", "rat_hgvsp"))
u <- merge(u, hum, by.x = c("human_symbol", "hgvs_p"), by.y = c("gene", "rat_hgvsp"), all.x = TRUE)

# --- CNV / LOH state per variant ---
cnv <- fread("/mnt/shared-workspace/shared/HeDe_cnv_segments.tsv")
setnames(cnv, c("ID", "chr", "start", "end", "nmark", "log2r", "width", "state"))
cnv[, `:=`(start = as.numeric(start), end = as.numeric(end), log2r = as.numeric(log2r))]
cnv[, chr := as.character(chr)]
get_cnv <- function(ch, p) {
  hit <- cnv[chr == ch & start <= p & end >= p]
  if (nrow(hit) == 0) return(list(NA_real_, NA_character_))
  list(hit$log2r[1], hit$state[1])
}
cnvinfo <- u[, { l <- get_cnv(chr, POS); .(cnv_log2r = l[[1]], cnv_state = l[[2]]) },
             by = seq_len(nrow(u))]
u[, `:=`(cnv_log2r = cnvinfo$cnv_log2r, cnv_state = cnvinfo$cnv_state)]
# LOH call: AF ~1 (homozygous) on a non-amplified segment
u[, loh := fifelse(AF >= 0.95 & cnv_state %in% c("neutral", "loss"),
                   fifelse(cnv_state == "loss", "homozygous_via_loss", "homozygous_CN-LOH"),
                   "heterozygous")]

# --- composite functional-evidence call ---
classify <- function(impact, sift, cadd, am, loh) {
  if (impact == "HIGH") return("high_impact_truncating")
  has_human <- !is.na(am) | !is.na(cadd)
  # human ortholog evidence takes precedence
  if (!is.na(am) & am == "likely_pathogenic") return("likely_damaging")
  if (!is.na(cadd) & cadd >= 25) return("likely_damaging")
  if (has_human & !is.na(am) & am == "likely_benign" & (is.na(cadd) | cadd < 20))
    return("likely_benign")
  if (has_human & !is.na(cadd) & cadd >= 20) return("possibly_damaging")
  # rat-SIFT-only evidence (no human score): lower confidence
  if (!is.na(sift) & grepl("deleterious", sift)) return("possibly_damaging")
  "uncertain"
}
u[, functional_call := mapply(classify, impact, rat_sift_pred, cadd_phred, am_class, loh)]

# order columns
setcolorder(u, c("human_symbol", "rat_symbol", "role", "is_hcc", "n_ctypes",
                 "CHROM", "POS", "REF", "ALT", "impact", "cons_group", "hgvs_p",
                 "GT", "AF", "rat_sift_pred", "rat_sift_score",
                 "human_hgvsp", "cadd_phred", "am_class", "am_score", "human_sift",
                 "cnv_state", "cnv_log2r", "loh", "functional_call"))
setorder(u, -is_hcc, human_symbol, POS)

fwrite(u, file.path(WD, "variants/cancer_mutations_functional_impact.tsv"), sep = "\t")

cat("=== rows:", nrow(u), "===\n")
cat("\nfunctional_call counts:\n"); print(u[, .N, by = functional_call][order(-N)])
cat("\nloh counts:\n"); print(u[, .N, by = loh][order(-N)])
cat("\n=== driver + deep-dive missense (key results) ===\n")
key <- u[human_symbol %in% c("TP53","KMT2D","EP300","FAT1","DICER1","ARID1B","POLQ","CLTC","TET1","KMT2C")]
print(key[, .(human_symbol, hgvs_p, AF, rat_sift_pred, human_hgvsp, cadd_phred, am_class, cnv_state, loh, functional_call)], width = 220)
