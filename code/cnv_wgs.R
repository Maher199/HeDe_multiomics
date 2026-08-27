#!/usr/bin/env Rscript
# Phase 5: WGS tumor/normal copy-number analysis, mirroring the PacBio method.
# mosdepth 20 kb bins -> log2(tumor/liver), each sample scaled to its chromosome-set
# median depth -> DNAcopy CBS (defaults) -> segment states at log2r +/-0.4 (the rule
# recovered from the PacBio segments) -> per-chromosome comparison with the PacBio
# karyotype table.
#
# Usage: Rscript cnv_wgs.R <tumor.20kb.regions.bed.gz> <normal.20kb.regions.bed.gz> <outdir>
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages({library(DNAcopy); library(readr); library(dplyr)})

args <- commandArgs(trailingOnly = TRUE)
tum_f <- args[1]; nor_f <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

CHRLEN <- "/mnt/results/deposition/zenodo/HeDe_PacBio_GRCr8_dataset/tables/GRCr8_chromosome_lengths.tsv"
PACKAR <- "/mnt/results/deposition/zenodo/HeDe_PacBio_GRCr8_dataset/tables/HeDe_karyotype_by_chromosome.tsv"

lens <- read_tsv(CHRLEN, show_col_types = FALSE) %>% rename(len = length)  # chr, accession, len
acc2chr <- setNames(lens$chr, lens$accession)

read_mos <- function(f) {
  d <- read_tsv(f, col_names = c("acc", "start", "end", "depth"),
                col_types = "ciid", progress = FALSE)
  d %>% filter(acc %in% names(acc2chr)) %>%
    mutate(chr = acc2chr[acc]) %>%
    select(chr, start, end, depth)
}

tum <- read_mos(tum_f); nor <- read_mos(nor_f)
cat("bins: tumor", nrow(tum), " normal", nrow(nor), "\n")

m <- inner_join(tum, nor, by = c("chr", "start", "end"), suffix = c("_tum", "_nor"))

med_t <- median(m$depth_tum); med_n <- median(m$depth_nor)
cat("median bin depth: tumor", med_t, " normal", med_n, "\n")

# reliability filter: normal bin depth >= 0.25x median (drops mappability deserts)
keep <- m$depth_nor >= 0.25 * med_n
cat("bins kept:", sum(keep), " dropped:", sum(!keep),
    sprintf("(%.2f%% dropped)\n", 100 * mean(!keep)))
mf <- m[keep, ]
mf$log2r <- log2((mf$depth_tum / med_t) / (mf$depth_nor / med_n))

chr_order <- lens$chr
mf$chr <- factor(mf$chr, levels = chr_order)
mf <- mf %>% arrange(chr, start)

cna <- CNA(genomdat = mf$log2r, chrom = as.character(mf$chr), maploc = mf$start,
           data.type = "logratio", sampleid = "HeDe_WGS_TvsN")
sm <- smooth.CNA(cna)
seg <- segment(sm, verbose = 0)
segs <- seg$output
segs$state <- ifelse(segs$seg.mean > 0.4, "gain", ifelse(segs$seg.mean < -0.4, "loss", "neutral"))
cat("segments:", nrow(segs), " states:",
    paste(names(table(segs$state)), table(segs$state), collapse = " / "), "\n")

seg_out <- segs %>%
  transmute(ID = "HeDe_WGS", chr = as.character(chrom), start = loc.start, end = loc.end,
            nmark = num.mark, log2r = round(seg.mean, 4),
            width = loc.end - loc.start, state)
write_tsv(seg_out, file.path(outdir, "HeDe_WGS_cnv_segments.tsv"))
write.csv(seg_out, file.path(outdir, "HeDe_WGS_cnv_segments.csv"), row.names = FALSE)

per_chr <- seg_out %>%
  group_by(chr) %>%
  summarise(
    n_seg = n(),
    width_gain = sum(width[state == "gain"]),
    width_loss = sum(width[state == "loss"]),
    width_neut = sum(width[state == "neutral"]),
    wlog2r = round(weighted.mean(log2r, width), 3),
    .groups = "drop") %>%
  left_join(lens, by = "chr") %>%
  mutate(covered = width_gain + width_loss + width_neut,
         frac_gain = round(width_gain / len, 3),
         frac_loss = round(width_loss / len, 3),
         frac_neut = round(width_neut / len, 3),
         call = case_when(n_seg == 0 ~ "absent",
                          frac_gain >= 0.50 ~ "gain",
                          frac_loss >= 0.50 ~ "loss",
                          frac_gain >= 0.15 ~ "partial gain",
                          frac_loss >= 0.15 ~ "partial loss",
                          TRUE ~ "neutral"))
write_tsv(per_chr, file.path(outdir, "HeDe_WGS_karyotype_by_chromosome.tsv"))

pac <- read_tsv(PACKAR, show_col_types = FALSE)
cmp <- pac %>%
  select(chr, n_seg, wlog2r, frac_gain, frac_loss, call) %>%
  rename_with(~ paste0(.x, "_pacbio"), .cols = -chr) %>%
  left_join(per_chr %>%
              select(chr, n_seg, wlog2r, frac_gain, frac_loss, call) %>%
              rename_with(~ paste0(.x, "_wgs"), .cols = -chr),
            by = "chr") %>%
  mutate(call_concordant = call_pacbio == call_wgs,
         dwlog2r = round(wlog2r_wgs - wlog2r_pacbio, 3))
write_tsv(cmp, file.path(outdir, "karyotype_pacbio_vs_wgs.tsv"))

cat("\n=== per-chromosome call comparison (PacBio vs WGS) ===\n")
print(as.data.frame(cmp[, c("chr", "call_pacbio", "call_wgs", "call_concordant",
                            "wlog2r_pacbio", "wlog2r_wgs", "dwlog2r")]), row.names = FALSE)
cat("\nconcordant calls:", sum(cmp$call_concordant), "/", nrow(cmp), "\n")

pc <- list(
  chr13_gain  = cmp$call_wgs[cmp$chr == "chr13"] == "gain",
  chr12_loss  = cmp$call_wgs[cmp$chr == "chr12"] == "loss",
  chrY_loss   = cmp$call_wgs[cmp$chr == "chrY"] == "loss",
  chr5_ploss  = cmp$call_wgs[cmp$chr == "chr5"] == "partial loss"
)
cat("\n=== positive controls ===\n")
for (k in names(pc)) cat(sprintf("  %-12s %s\n", k, ifelse(isTRUE(pc[[k]]), "PASS", "FAIL/CHECK")))

chr5 <- seg_out %>% filter(chr == "chr5") %>% arrange(log2r)
cat("\ndeepest chr5 segment (Cdkn2a/b region expected):\n")
print(as.data.frame(head(chr5, 3)), row.names = FALSE)
cat("\nDONE\n")
