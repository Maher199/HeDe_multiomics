#!/usr/bin/env Rscript
# Scientific Data Data Descriptor figures for the HeDe PacBio HiFi dataset.
# Multi-panel, lower-case bold a/b/c labels, Liberation Sans, white bg, Okabe-Ito palette.
suppressMessages({library(data.table); library(ggplot2); library(svglite); library(patchwork)})

WD  <- "/workspace/hede_followup"
FT  <- "/mnt/results/followup_tables"
DEP <- "/mnt/results/deposition/zenodo/HeDe_PacBio_GRCr8_dataset"
COV <- "/mnt/shared-workspace/shared/coverage"
OUT <- "/workspace/hede_followup/dd_figures"; dir.create(OUT, showWarnings = FALSE)
FONT <- "Liberation Sans"
OI <- c(grey="#999999", orange="#E69F00", sky="#56B4E9", green="#009E73",
        yellow="#F0E442", blue="#0072B2", verm="#D55E00", pink="#CC79A7", black="#000000")

theme_pub <- theme_classic(base_size = 9, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.3),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold", hjust = -0.05),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(5, 7, 5, 5))

saveboth <- function(p, name, w, h) {
  svglite(file.path(OUT, paste0(name, ".svg")), width = w, height = h, fix_text_size = FALSE)
  print(p); dev.off()
  png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
}

# NC_ accession -> chromosome label (GRCr8): NC_0860XX -> chr(XX-18); 39->X, 40->Y
nc2chr <- function(acc) {
  n <- as.integer(sub("NC_0860([0-9]+)\\.1", "\\1", acc))
  ifelse(n >= 19 & n <= 38, paste0("chr", n - 18),
         ifelse(n == 39, "chrX", ifelse(n == 40, "chrY", NA)))
}

## ============================================================
## FIG 1: Sequencing & alignment QC
##  a) HiFi read-length distribution (N50 marked)
##  b) mean coverage per chromosome
##  c) cumulative genome coverage
## ============================================================
# Read lengths for ALL 5,565,830 HiFi reads. Regenerate from the raw BAM index with
#   Rscript code/read_lengths_from_pbi.R H_mix.hifi_reads.bam.pbi tables/
# Do not substitute a head-of-file subsample here: the first reads in BAM order are
# systematically longer (early ZMWs), which inflates the N50 (first 500,000 reads give
# 23,526 bp versus 23,294 bp for the full set).
rl <- fread(file.path(DEP, "tables/read_lengths_full.txt.gz"), header = FALSE, col.names = "len")
n50 <- {
  l <- sort(as.numeric(rl$len), decreasing = TRUE); cs <- cumsum(l); l[cs >= sum(l)/2][1]
}
maxc <- max(hist(rl$len / 1000, breaks = seq(0, ceiling(max(rl$len / 1000)) + 1, by = 0.5), plot = FALSE)$counts)
p1a <- ggplot(rl, aes(x = len / 1000)) +
  geom_histogram(binwidth = 0.5, fill = OI[["blue"]], colour = NA) +
  geom_vline(xintercept = n50 / 1000, linetype = "dashed", linewidth = 0.7, colour = OI[["verm"]]) +
  annotate("text", x = (n50 / 1000) + 1.2, y = maxc * 1.08,
           label = paste0("N50 = ", round(n50/1000, 1), " kb"),
           hjust = 0, vjust = 1, size = 3, family = FONT, fontface = "bold", colour = OI[["verm"]]) +
  coord_cartesian(xlim = c(0, 45)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "HiFi read length (kb)", y = "Reads", title = "a") + theme_pub

ms <- fread(file.path(COV, "HeDe_500bp.mosdepth.summary.txt"))
ms <- ms[!grepl("_region", chrom)]
ms[, chr := nc2chr(chrom)]
ms <- ms[!is.na(chr)]
ms[, chr := factor(chr, levels = c(paste0("chr", 1:20), "chrX", "chrY"))]
genome_mean <- weighted.mean(ms$mean, ms$length)
p1b <- ggplot(ms, aes(x = chr, y = mean)) +
  geom_col(width = 0.72, fill = OI[["blue"]]) +
  geom_hline(yintercept = genome_mean, linetype = "dashed", linewidth = 0.35, colour = OI[["verm"]]) +
  annotate("text", x = 0.4, y = Inf, label = paste0("mean = ", round(genome_mean, 1), "x"),
           hjust = 0, vjust = 1.5, size = 2.7, family = FONT, fontface = "bold", colour = OI[["verm"]]) +
  labs(x = "Chromosome", y = "Mean coverage (x)", title = "b") +
  theme_pub + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5))

gd <- fread(file.path(COV, "HeDe_500bp.mosdepth.global.dist.txt"), header = FALSE,
            col.names = c("chrom", "cov", "prop"))
gd <- gd[chrom == "total"][order(cov)]
p1c <- ggplot(gd, aes(x = cov, y = prop)) +
  geom_line(linewidth = 0.6, colour = OI[["green"]]) +
  geom_vline(xintercept = round(genome_mean), linetype = "dashed", linewidth = 0.35, colour = "grey35") +
  coord_cartesian(xlim = c(0, 120)) +
  labs(x = "Coverage threshold (x)", y = "Fraction of genome ≥ threshold", title = "c") +
  theme_pub

fig1 <- p1a + p1b + p1c + plot_layout(widths = c(1.15, 1.15, 1))
saveboth(fig1, "dd_fig1_sequencing_qc", 7.4, 2.7)

## ============================================================
## FIG 2: De novo assembly quality
##  a) cumulative contig length (Nx curve, N50 marked)
##  b) contig-length distribution (top contigs)
##  c) BUSCO completeness
## ============================================================
ct <- fread(file.path(WD, "dd_figdata/contig_lengths.tsv"))
ct <- ct[order(-length)]
ct[, rank := .I]
ct[, cum := cumsum(as.numeric(length))]
tot <- sum(ct$length)
ct[, frac := cum / tot]
n50c <- ct[frac >= 0.5][1]
n50rank <- which(ct$frac >= 0.5)[1]
p2a <- ggplot(ct, aes(x = rank, y = frac)) +
  geom_line(linewidth = 0.6, colour = OI[["blue"]]) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_point(data = ct[n50rank], aes(x = rank, y = frac), colour = OI[["verm"]], size = 2) +
  annotate("text", x = n50rank, y = 0.43,
           label = paste0("N50 = ", round(n50c$length/1e6, 1), " Mb\n(contig ", n50rank, ")"),
           hjust = 0.5, vjust = 1, size = 2.6, family = FONT, colour = OI[["verm"]]) +
  scale_x_log10(expand = expansion(mult = c(0.01, 0.05))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Contig rank (log10, largest first)", y = "Cumulative fraction of assembly", title = "a") +
  theme_pub

p2b <- ggplot(ct, aes(x = length)) +
  geom_histogram(bins = 50, fill = OI[["sky"]], colour = NA) +
  scale_x_log10(labels = function(x) paste0(format(x/1e6, trim = TRUE), " Mb")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Contig length (log10)", y = "Contigs", title = "b") +
  theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

busco <- data.table(
  cat = factor(c("Complete single-copy", "Complete duplicated", "Fragmented", "Missing"),
               levels = c("Complete single-copy", "Complete duplicated", "Fragmented", "Missing")),
  pct = c(98.28, 1.51, 0.08, 0.14))
p2c <- ggplot(busco, aes(x = "", y = pct, fill = cat)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pct > 2, paste0(round(pct, 1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 2.6, family = FONT, colour = "white") +
  scale_fill_manual(values = c("Complete single-copy" = OI[["blue"]],
                               "Complete duplicated" = OI[["sky"]],
                               "Fragmented" = OI[["orange"]],
                               "Missing" = OI[["verm"]])) +
  labs(x = NULL, y = "BUSCO genes (%)", fill = NULL, title = "c") +
  theme_pub + theme(legend.position = "right")

fig2 <- p2a + p2b + p2c + plot_layout(widths = c(1.2, 0.95, 0.7))
saveboth(fig2, "dd_fig2_assembly_quality", 7.4, 2.9)

## ============================================================
## FIG 3: Variant & structural-variant landscape
##  a) consequence spectrum (log)
##  b) SV counts by type
##  c) SV size distribution
##  d) 6-class mutational spectrum
## ============================================================
cons <- fread(file.path(FT, "consequence_summary.tsv"))
cons <- cons[, .(N = sum(N)), by = .(cons_group, impact)]
cons[, cons_group := factor(cons_group, levels = unique(cons[order(-N)]$cons_group))]
p3a <- ggplot(cons, aes(x = reorder(cons_group, N), y = N, fill = impact)) +
  geom_col(width = 0.72) + coord_flip() +
  scale_y_log10(expand = expansion(mult = c(0, 0.12)), labels = scales::comma) +
  scale_fill_manual(values = c(HIGH = OI[["verm"]], MODERATE = OI[["orange"]],
                               LOW = OI[["sky"]], MODIFIER = OI[["grey"]])) +
  labs(x = NULL, y = "Variants (log10)", fill = "Impact", title = "a") +
  theme_pub + theme(legend.position = c(0.72, 0.32))

sv <- fread(file.path(DEP, "sv/HeDe_sv_table.tsv"), header = FALSE,
            col.names = c("chrom", "pos", "type", "len"))
svN <- sv[, .N, by = type][order(-N)]
svN[, type := factor(type, levels = svN$type)]
p3b <- ggplot(svN, aes(x = reorder(type, N), y = N, fill = type)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = format(N, big.mark = ",")), hjust = -0.1, size = 2.6, family = FONT) +
  coord_flip() +
  scale_y_log10(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = c(DEL = OI[["verm"]], INS = OI[["blue"]], DUP = OI[["green"]],
                               INV = OI[["orange"]], BND = OI[["grey"]], INVDUP = OI[["pink"]])) +
  labs(x = NULL, y = "SVs (log10)", title = "b") + theme_pub

svs <- sv[abs(len) >= 50]
svs[, size := abs(len)]
p3c <- ggplot(svs, aes(x = size)) +
  geom_histogram(bins = 60, fill = OI[["orange"]], colour = NA) +
  scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "SV size (bp, log10)", y = "SVs", title = "c") + theme_pub

sp <- fread(file.path(DEP, "tables/mutational_spectrum_6class.csv"))
sp[, subtype := factor(subtype, levels = c("C>A","C>G","C>T","T>A","T>C","T>G"))]
p3d <- ggplot(sp, aes(x = subtype, y = count, fill = subtype)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 2.3, family = FONT, angle = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), labels = function(x) format(x, big.mark=",", scientific=FALSE)) +
  scale_fill_manual(values = c("C>A"=OI[["sky"]], "C>G"=OI[["black"]], "C>T"=OI[["verm"]],
                               "T>A"=OI[["grey"]], "T>C"=OI[["green"]], "T>G"=OI[["pink"]])) +
  labs(x = "Substitution class", y = "Variants", title = "d") + theme_pub

fig3 <- (p3a + p3b) / (p3c + p3d)
saveboth(fig3, "dd_fig3_variant_sv_landscape", 7.4, 5.4)

## ============================================================
## FIG 5: Methylation landscape
##  a) CGI methylation states
##  b) retrotransposon family delta vs genome
##  c) imprinted-locus methylation
## ============================================================
cgi <- fread(file.path(FT, "cgi_5mC.tsv"))
cgiN <- cgi[, .N, by = meth_state]
cgiN[, meth_state := factor(meth_state, levels = c("hypomethylated", "intermediate", "hypermethylated"))]
p5a <- ggplot(cgiN, aes(x = meth_state, y = N, fill = meth_state)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = format(N, big.mark=",")), vjust = -0.4, size = 2.6, family = FONT) +
  scale_fill_manual(values = c(hypomethylated = OI[["sky"]], intermediate = OI[["grey"]],
                               hypermethylated = OI[["verm"]])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "CpG-island state", y = "Islands", title = "a") +
  theme_pub + theme(axis.text.x = element_text(angle = 20, hjust = 1))

rep <- fread(file.path(FT, "retrotransposon_5mC.tsv"))
repM <- rep[n_cpg >= 10000]
repM[, fam := sub("^(LINE|SINE|LTR)/", "", classfam)]
repM[, fam := paste0(fam, " (", rclass, ")")]
p5b <- ggplot(repM, aes(x = reorder(fam, delta_vs_genome), y = delta_vs_genome, fill = rclass)) +
  geom_col(width = 0.7) + coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_y_continuous(breaks = c(-0.05, 0, 0.05, 0.10), limits = c(-0.06, 0.135)) +
  scale_fill_manual(values = c(LINE = OI[["blue"]], SINE = OI[["orange"]], LTR = OI[["green"]])) +
  labs(x = NULL, y = expression(Delta*"5mC vs genome"), fill = "Class", title = "b") +
  theme_pub + theme(legend.position = c(0.82, 0.28))

imp <- fread(file.path(FT, "imprinted_5mC.tsv"))
imp[, locus := sub("_promDMR$", "", region_id)]
imp[, locus := sub("_(ICR|IGDMR)$", " \\(\\1\\)", locus)]
imp[, loi_state := factor(loi_state, levels = c("loss(hypo)", "intact(~0.5)", "gain(hyper)"))]
p5c <- ggplot(imp, aes(x = reorder(locus, mean5mC), y = mean5mC, fill = loi_state)) +
  geom_point(shape = 21, size = 2.4, colour = "black", stroke = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  coord_flip() + scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = c("loss(hypo)" = OI[["sky"]], "intact(~0.5)" = OI[["grey"]],
                               "gain(hyper)" = OI[["verm"]])) +
  labs(x = NULL, y = "Mean 5mC", fill = "Imprinting", title = "c") +
  theme_pub + theme(axis.text.y = element_text(size = 7),
                    legend.position = "bottom", legend.direction = "horizontal")

fig5 <- (p5a + p5b) / p5c + plot_layout(heights = c(1, 1.4))
saveboth(fig5, "dd_fig5_methylation_landscape", 7.4, 6.4)

## ============================================================
## FIG 6: Technical validation
##  a) het-SNP AF density (Tp53 vs chr10 vs genome) — no LOH
##  b) TP53 R271S read counts
##  c) Cdkn2a/Cdkn2b homozygous deletion coverage
## ============================================================
af <- readRDS(file.path(WD, "variants/tp53_loh_af_vectors.rds"))
tp <- af$tp53; c10 <- af$chr10_rest; gw <- af$genome_other
set.seed(1); if (length(gw) > 200000) gw <- sample(gw, 200000)
df <- rbind(data.table(AF = tp, region = "Tp53 locus"),
            data.table(AF = c10, region = "chr10 (rest)"),
            data.table(AF = gw, region = "Genome (other)"))
df[, region := factor(region, levels = c("Tp53 locus", "chr10 (rest)", "Genome (other)"))]
p6a <- ggplot(df, aes(x = AF, colour = region)) +
  geom_density(linewidth = 0.6, adjust = 1.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  scale_color_manual(values = c("Tp53 locus" = OI[["verm"]], "chr10 (rest)" = OI[["blue"]],
                                "Genome (other)" = OI[["grey"]])) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Heterozygous-SNP allele fraction", y = "Density", colour = NULL, title = "a") +
  theme_pub + theme(legend.position = c(0.18, 0.72))

r271 <- data.table(allele = c("C (ref)", "A (R271S)"), reads = c(20, 28))
p6b <- ggplot(r271, aes(x = allele, y = reads, fill = allele)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = reads), vjust = -0.4, size = 3, family = FONT) +
  scale_fill_manual(values = c("C (ref)" = OI[["grey"]], "A (R271S)" = OI[["verm"]])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Read count", title = "b") + theme_pub

cov <- fread(file.path(WD, "variants/cdkn2a_coverage_500bp.tsv"), col.names = c("pos", "cov"))
bp <- fread(file.path(WD, "variants/cdkn2a_breakpoint_confirmation.tsv"))
bpv <- function(k) bp[feature == k]$value
bl <- as.numeric(bpv("breakpoint_left")); br <- as.numeric(bpv("breakpoint_right"))
genes <- data.table(gene = c("Cdkn2a", "Cdkn2b"),
                    start = c(109100763, 109123308), end = c(109114448, 109134906))
flank_mean <- as.numeric(bpv("coverage_flank_mean"))
p6c <- ggplot(cov, aes(x = pos / 1e6, y = cov)) +
  annotate("rect", xmin = bl/1e6, xmax = br/1e6, ymin = -Inf, ymax = Inf,
           fill = OI[["verm"]], alpha = 0.12) +
  geom_line(linewidth = 0.35, colour = OI[["blue"]]) +
  geom_hline(yintercept = flank_mean, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_vline(xintercept = c(bl, br)/1e6, linetype = "dotted", linewidth = 0.4, colour = OI[["verm"]]) +
  geom_segment(data = genes, aes(x = start/1e6, xend = end/1e6, y = -4, yend = -4),
               linewidth = 3, colour = OI[["green"]], inherit.aes = FALSE) +
  geom_text(data = genes, aes(x = (start+end)/2/1e6, y = -7.5, label = gene),
            size = 2.5, family = FONT, fontface = "italic", inherit.aes = FALSE) +
  coord_cartesian(ylim = c(-9, max(cov$cov)*1.02), clip = "off") +
  labs(x = "chr5 position (Mb)", y = "HiFi read depth", title = "c") + theme_pub

fig6 <- (p6a + p6b) / p6c + plot_layout(heights = c(1, 1))
saveboth(fig6, "dd_fig6_technical_validation", 7.4, 5.4)

cat("DD figures written to", OUT, "\n")
print(list.files(OUT, pattern = "\\.svg$"))
