#!/usr/bin/env Rscript
# Analysis C1: Tp53 LOH via heterozygous-SNP allele-fraction (AF) distribution.
# Logic: in a diploid region without LOH, germline het SNPs have AF ~ 0.5.
#   - copy-neutral LOH -> deficit of het SNPs (runs of homozygosity) + AF shifts to 0/1
#   - deletion LOH     -> loss of het SNPs + copy-number loss (Tp53 is diploid, so ruled out)
# Compare Tp53 locus AF distribution & het density vs rest of chr and genome-wide baseline.
suppressMessages({ library(data.table) })
WD <- "/workspace/hede_followup/variants"

gw <- fread(file.path(WD, "gw_het_af.tsv"), sep = "\t", header = FALSE,
            col.names = c("chrom", "pos", "af"))
# drop multi-allelic (comma-separated AF) and non-numeric; keep simple het sites
gw <- gw[!grepl(",", af)]
gw[, af := as.numeric(af)]
gw <- gw[!is.na(af) & af > 0 & af < 1]
message("genome-wide het sites (biallelic, 0<AF<1): ", nrow(gw))

TP53_CHR <- "NC_086028.1"
TP53_G0 <- 54798871; TP53_G1 <- 54810300
WIN <- 2e6
L0 <- TP53_G0 - WIN; L1 <- TP53_G1 + WIN   # 52.8 - 56.8 Mb locus window

gw[, region := fcase(
  chrom == TP53_CHR & pos >= L0 & pos <= L1, "Tp53_locus",
  chrom == TP53_CHR,                          "chr10_rest",
  default = "genome_other")]
reg_tab <- gw[, .(
  n_het = .N,
  median_AF = round(median(af), 4),
  mean_AF = round(mean(af), 4),
  sd_AF = round(sd(af), 4),
  frac_AF_0.4_0.6 = round(mean(af >= 0.4 & af <= 0.6), 4),
  frac_AF_extreme = round(mean(af < 0.25 | af > 0.75), 4)   # LOH shoulder signal
), by = region]
# het density per Mb
chr_len <- gw[, .(span = max(pos)), by = chrom]
tp53_span_mb <- (L1 - L0) / 1e6
chr10_rest_mb <- (chr_len[chrom == TP53_CHR]$span - (L1 - L0)) / 1e6
other_mb <- (sum(chr_len$span) - chr_len[chrom == TP53_CHR]$span) / 1e6
reg_tab[, het_per_Mb := round(fcase(
  region == "Tp53_locus", n_het / tp53_span_mb,
  region == "chr10_rest",  n_het / chr10_rest_mb,
  region == "genome_other", n_het / other_mb), 1)]
reg_tab[, ord := match(region, c("Tp53_locus", "chr10_rest", "genome_other"))]
setorder(reg_tab, ord); reg_tab[, ord := NULL]
print(reg_tab)
fwrite(reg_tab, file.path(WD, "tp53_loh_af_summary.tsv"), sep = "\t")

# --- tests ---
tp <- gw[region == "Tp53_locus"]$af
other <- gw[region == "genome_other"]$af
chr8r <- gw[region == "chr10_rest"]$af
cat("\n-- one-sample Wilcoxon: Tp53 locus AF vs 0.5 --\n")
print(wilcox.test(tp, mu = 0.5))
cat("\n-- KS test: Tp53 locus AF vs genome-wide baseline --\n")
print(ks.test(tp, other))
cat("\n-- KS test: Tp53 locus AF vs rest of chr10 --\n")
print(ks.test(tp, chr8r))

# --- R273S zygosity: AF=0.5833 (28/48). Binomial test vs 0.5 ---
cat("\n-- R273S (C>A, AD 19/28) binomial test vs AF=0.5 --\n")
print(binom.test(28, 48, p = 0.5))

# save AF vectors for plotting later
saveRDS(list(tp53 = tp, chr10_rest = chr8r, genome_other = other),
        file.path(WD, "tp53_loh_af_vectors.rds"))
cat("\nDONE\n")
