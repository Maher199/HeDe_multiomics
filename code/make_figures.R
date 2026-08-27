#!/usr/bin/env Rscript
# Deliverable D: publication-ready multi-panel SVG figures for the HeDe PacBio rat-HCC study.
# ggplot2 + svglite, Liberation Sans, colorblind-safe (Okabe-Ito) palette.
suppressMessages({library(data.table); library(ggplot2); library(svglite); library(patchwork)})

WD  <- "/workspace/hede_followup"
OUT <- "/workspace/hede_followup/figures"; dir.create(OUT, showWarnings = FALSE)
FONT <- "Liberation Sans"
# Okabe-Ito colorblind-safe palette
OI <- c(grey="#999999", orange="#E69F00", sky="#56B4E9", green="#009E73",
        yellow="#F0E442", blue="#0072B2", verm="#D55E00", pink="#CC79A7", black="#000000")
HEX <- list(grey="#999999", orange="#E69F00", sky="#56B4E9", green="#009E73",
            yellow="#F0E442", blue="#0072B2", verm="#D55E00", pink="#CC79A7", black="#000000")

theme_pub <- theme_classic(base_size = 9, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.3),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        plot.margin = margin(4, 6, 4, 4))

savefig <- function(p, name, w, h) {
  svglite(file.path(OUT, paste0(name, ".svg")), width = w, height = h, fix_text_size = FALSE)
  print(p); dev.off()
}

## ============================================================
## FIGURE 1: Coding-mutation landscape
##  A) consequence spectrum (primary annotation, log scale)
##  B) high-impact cancer-gene mutations (unique variants)
## ============================================================
cons <- fread(file.path(WD, "variants/consequence_summary.tsv"))
# aggregate by cons_group x impact (avoid duplicate-level factor issues)
cons <- cons[, .(N = sum(N)), by = .(cons_group, impact)]
cons[, cons_group := factor(cons_group, levels = unique(cons[order(-N)]$cons_group))]
pA <- ggplot(cons, aes(x = reorder(cons_group, N), y = N, fill = impact)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_y_log10(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = c(HIGH = HEX$verm, MODERATE = HEX$orange,
                               LOW = HEX$sky, MODIFIER = HEX$grey)) +
  labs(x = NULL, y = "Novel variants (log10)", fill = "Impact", title = "A") +
  theme_pub + theme(legend.position = c(0.75, 0.35))

# high-impact cancer genes (unique variants, dedup by gene+pos)
cg <- fread(file.path(WD, "variants/cancer_coding_mutations.tsv"))
hi <- cg[impact == "HIGH"]
hiu <- unique(hi, by = c("human_symbol", "CHROM", "POS"))
# one row per unique variant; label gene + consequence
hiu[, glabel := factor(human_symbol, levels = hiu[order(AF)]$human_symbol)]
hiu[, cons_short := sub("/start-stop", "", cons_group)]
pB <- ggplot(hiu, aes(x = reorder(human_symbol, AF), y = AF, fill = cons_short)) +
  geom_point(shape = 21, size = 2.6, color = "black", stroke = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
  scale_fill_manual(values = c("nonsense" = HEX$verm, "frameshift" = HEX$blue,
                               "splice" = HEX$green)) +
  labs(x = NULL, y = "Variant allele fraction", fill = "Consequence",
       title = "B") +
  theme_pub + theme(legend.position = "bottom", legend.direction = "horizontal")

fig1 <- pA + pB + plot_layout(widths = c(1, 1))
savefig(fig1, "fig1_coding_mutation_landscape", 7.2, 3.4)

## ============================================================
## FIGURE 2: Tp53 locus — no LOH
##  A) het-SNP AF density: Tp53 locus vs chr10 rest vs genome
##  B) R271S read counts (balanced het)
## ============================================================
af <- readRDS(file.path(WD, "variants/tp53_loh_af_vectors.rds"))
tp <- af$tp53; c10 <- af$chr10_rest; gw <- af$genome_other
# subsample genome for plotting speed
set.seed(1); if (length(gw) > 200000) gw <- sample(gw, 200000)
df <- rbind(data.table(AF = tp, region = "Tp53 locus"),
            data.table(AF = c10, region = "chr10 (rest)"),
            data.table(AF = gw, region = "Genome (other)"))
df[, region := factor(region, levels = c("Tp53 locus", "chr10 (rest)", "Genome (other)"))]
pA2 <- ggplot(df, aes(x = AF, color = region)) +
  geom_density(linewidth = 0.6, adjust = 1.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  scale_color_manual(values = c("Tp53 locus" = HEX$verm, "chr10 (rest)" = HEX$blue,
                                "Genome (other)" = HEX$grey)) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Heterozygous-SNP allele fraction", y = "Density", color = NULL, title = "A") +
  theme_pub + theme(legend.position = c(0.72, 0.7))

# R271S read counts (from breakpoint confirmation: 20 C / 28 A)
r271 <- data.table(allele = c("C (ref)", "A (R271S)"), reads = c(20, 28))
pB2 <- ggplot(r271, aes(x = allele, y = reads, fill = allele)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = reads), vjust = -0.4, size = 3, family = FONT) +
  scale_fill_manual(values = c("C (ref)" = HEX$grey, "A (R271S)" = HEX$verm)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Read count", title = "B") +
  theme_pub

fig2 <- pA2 + pB2 + plot_layout(widths = c(1.5, 1))
savefig(fig2, "fig2_tp53_no_loh", 7.2, 2.9)

## ============================================================
## FIGURE 3: Cdkn2a/Cdkn2b (INK4/ARF) homozygous deletion
##  coverage profile + gene models + breakpoints
## ============================================================
cov <- fread(file.path(WD, "variants/cdkn2a_coverage_500bp.tsv"),
             col.names = c("pos", "cov"))
bp <- fread(file.path(WD, "variants/cdkn2a_breakpoint_confirmation.tsv"))
bpv <- function(k) bp[feature == k]$value
bl <- as.numeric(bpv("breakpoint_left")); br <- as.numeric(bpv("breakpoint_right"))
# gene models (GRCr8 coords)
genes <- data.table(
  gene = c("Cdkn2a", "Cdkn2b"),
  start = c(109100763, 109123308), end = c(109114448, 109134906), y = c(1, 1))
del <- data.table(xmin = bl, xmax = br)
flank_mean <- as.numeric(bpv("coverage_flank_mean"))

p3 <- ggplot(cov, aes(x = pos / 1e6, y = cov)) +
  annotate("rect", xmin = bl/1e6, xmax = br/1e6, ymin = -Inf, ymax = Inf,
           fill = HEX$verm, alpha = 0.12) +
  geom_line(linewidth = 0.35, color = HEX$blue) +
  geom_hline(yintercept = flank_mean, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_vline(xintercept = c(bl, br)/1e6, linetype = "dotted", linewidth = 0.4, color = HEX$verm) +
  geom_segment(data = genes, aes(x = start/1e6, xend = end/1e6, y = -4, yend = -4),
               linewidth = 3, color = HEX$green, inherit.aes = FALSE) +
  geom_text(data = genes, aes(x = (start+end)/2/1e6, y = -7, label = gene),
            size = 2.6, family = FONT, fontface = "italic", inherit.aes = FALSE) +
  coord_cartesian(ylim = c(-9, max(cov$cov)*1.02), clip = "off") +
  labs(x = "chr5 position (Mb)", y = "HiFi read depth", title = NULL) +
  theme_pub
savefig(p3, "fig3_cdkn2a_deletion", 7.2, 2.7)

## ============================================================
## FIGURE 4: Methylation landscape
##  A) CGI methylation states (count by state)
##  B) retrotransposon family methylation (delta vs genome)
##  C) imprinted loci LOI (mean5mC by locus)
## ============================================================
cgi <- fread(file.path(WD, "methyl/cgi_5mC.tsv"))
cgiN <- cgi[, .N, by = meth_state]
cgiN[, meth_state := factor(meth_state, levels = c("hypomethylated", "intermediate", "hypermethylated"))]
pA4 <- ggplot(cgiN, aes(x = meth_state, y = N, fill = meth_state)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = N), vjust = -0.4, size = 2.7, family = FONT) +
  scale_fill_manual(values = c(hypomethylated = HEX$sky, intermediate = HEX$grey,
                               hypermethylated = HEX$verm)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "CpG-island state", y = "No. islands", title = "A") +
  theme_pub + theme(axis.text.x = element_text(angle = 20, hjust = 1))

rep <- fread(file.path(WD, "methyl/retrotransposon_5mC.tsv"))
repM <- rep[n_cpg >= 10000]  # well-covered families
repM[, fam := sub("^(LINE|SINE|LTR)/", "", classfam)]
repM[, fam := paste0(fam, " (", rclass, ")")]
pB4 <- ggplot(repM, aes(x = reorder(fam, delta_vs_genome), y = delta_vs_genome, fill = rclass)) +
  geom_col(width = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(LINE = HEX$blue, SINE = HEX$orange, LTR = HEX$green)) +
  labs(x = NULL, y = expression(Delta*"5mC vs genome"), fill = "Class", title = "B") +
  theme_pub + theme(legend.position = c(0.8, 0.25))

imp <- fread(file.path(WD, "methyl/imprinted_5mC.tsv"))
# disambiguate locus names (keep DMR type for duplicated genes like Dlk1)
imp[, locus := sub("_promDMR$", "", region_id)]
imp[, locus := sub("_(ICR|IGDMR)$", " \\(\\1\\)", locus)]
imp[, loi_state := factor(loi_state, levels = c("loss(hypo)", "intact(~0.5)", "gain(hyper)"))]
pC4 <- ggplot(imp, aes(x = reorder(locus, mean5mC), y = mean5mC, fill = loi_state)) +
  geom_point(shape = 21, size = 2.4, color = "black", stroke = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = c("loss(hypo)" = HEX$sky, "intact(~0.5)" = HEX$grey,
                               "gain(hyper)" = HEX$verm)) +
  labs(x = NULL, y = "Mean 5mC", fill = "Imprinting", title = "C") +
  theme_pub + theme(axis.text.y = element_text(size = 7),
                    legend.position = "bottom", legend.direction = "horizontal")

fig4 <- (pA4 + pB4) / pC4 + plot_layout(heights = c(1, 1.5))
savefig(fig4, "fig4_methylation_landscape", 7.2, 6.8)

cat("Figures written to ", OUT, "\n")
print(list.files(OUT, pattern = "\\.svg$"))
