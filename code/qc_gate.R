#!/usr/bin/env Rscript
# Phase 2 QC gate: parse flagstat + mosdepth summaries for both WGS samples,
# check acceptance criteria (>=90% mapped, >=25x mean coverage), infer the liver
# donor's sex from chrX/chrY coverage (needed for chrY-loss interpretation),
# and write a QC table for the report.
# Usage: Rscript qc_gate.R <bam_dir_root> <outdir>
#   expects <bam_dir_root>/<SAMPLE>/<SAMPLE>.flagstat.txt and .20kb.mosdepth.summary.txt
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages({library(readr); library(dplyr); library(tidyr)})

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]; outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

CHRLEN <- "/mnt/results/deposition/zenodo/HeDe_PacBio_GRCr8_dataset/tables/GRCr8_chromosome_lengths.tsv"
lens <- read_tsv(CHRLEN, show_col_types = FALSE)
genome_bp <- sum(lens$len)

parse_flagstat <- function(f) {
  l <- readLines(f)
  get <- function(pat) as.numeric(sub("^ *([0-9]+) .*", "\\1", grep(pat, l, value = TRUE)[1]))
  tibble(total_reads = get("in total"),
         mapped = get("mapped \\("),
         proper_pairs = get("properly paired"),
         duplicates = get("duplicates"))
}

parse_mosdepth_summary <- function(f) {
  d <- read_tsv(f, show_col_types = FALSE)   # chrom, length, bases, mean, min, max
  d
}

samples <- c("H_1", "RL_1")
qc <- lapply(samples, function(s) {
  fs <- parse_flagstat(file.path(root, s, paste0(s, ".flagstat.txt")))
  ms <- parse_mosdepth_summary(file.path(root, s, paste0(s, ".20kb.mosdepth.summary.txt")))
  ms_chr <- ms %>% filter(chrom %in% lens$accession)
  tot <- ms %>% filter(chrom == "total")
  # weighted mean over placed chromosomes only
  wmean <- sum(ms_chr$mean * ms_chr$length) / sum(ms_chr$length)
  chrx <- ms_chr$mean[ms_chr$chrom == "NC_086039.1"]
  chry <- ms_chr$mean[ms_chr$chrom == "NC_086040.1"]
  auto <- mean(ms_chr$mean[!ms_chr$chrom %in% c("NC_086039.1", "NC_086040.1")])
  fs %>% mutate(sample = s,
                pct_mapped = round(100 * mapped / total_reads, 2),
                pct_dup = round(100 * duplicates / total_reads, 2),
                mean_cov_total = round(tot$mean[1], 2),
                mean_cov_chrs = round(wmean, 2),
                cov_autosomes = round(auto, 2),
                cov_chrX = round(chrx, 2),
                cov_chrY = round(chry, 2),
                chrY_over_auto = round(chry / auto, 3),
                chrX_over_auto = round(chrx / auto, 3))
}) %>% bind_rows()

qc <- qc %>% mutate(
  sex_inference = case_when(
    chrY_over_auto > 0.4 & chrX_over_auto < 0.75 ~ "male",
    chrY_over_auto < 0.05 & chrX_over_auto > 0.85 ~ "female",
    TRUE ~ "ambiguous"),
  pass_mapped = pct_mapped >= 90,
  pass_cov = mean_cov_chrs >= 25)

write_tsv(qc, file.path(outdir, "wgs_qc_gate.tsv"))
print(as.data.frame(qc), row.names = FALSE)

cat("\n=== gates ===\n")
cat("mapped >=90%:  ", paste(qc$sample, qc$pass_mapped), "\n")
cat("coverage >=25x:", paste(qc$sample, qc$pass_cov), "\n")
cat("liver (RL_1) sex:", qc$sex_inference[qc$sample == "RL_1"],
    " (chrY/auto =", qc$chrY_over_auto[qc$sample == "RL_1"], ")\n")
cat("tumor (H_1) sex:", qc$sex_inference[qc$sample == "H_1"],
    " (chrY/auto =", qc$chrY_over_auto[qc$sample == "H_1"], ")\n")
