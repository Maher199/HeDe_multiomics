#!/usr/bin/env Rscript
# Figure: loss-of-function vs gain-of-function cancer-gene mutations in HeDe.
# LoF = truncating/splice mutations in tumour suppressors; GoF = mutations in oncogenes.
suppressMessages({library(data.table); library(ggplot2); library(svglite); library(patchwork)})

WD  <- "/workspace/hede_followup"
OUT <- "/workspace/hede_followup/figures"; dir.create(OUT, showWarnings = FALSE)
FONT <- "Liberation Sans"
HEX <- list(grey="#999999", orange="#E69F00", sky="#56B4E9", green="#009E73",
            yellow="#F0E442", blue="#0072B2", verm="#D55E00", pink="#CC79A7", black="#000000")
theme_pub <- theme_classic(base_size = 9, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.3),
        axis.title = element_text(size = 9), axis.text = element_text(size = 8),
        plot.title = element_text(size = 9.5, face = "bold"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"), plot.margin = margin(4,6,4,4))

cc <- fread(file.path(WD, "variants/cancer_coding_mutations.tsv"))
# dedup to unique variants (gene + pos)
u <- unique(cc, by = c("human_symbol","CHROM","POS"))
u[, cons_short := sub("/start-stop", "", cons_group)]
# functional consequence class for shaping
u[, func := fifelse(impact == "HIGH", "truncating/splice",
            fifelse(cons_group == "missense", "missense", "in-frame indel"))]
u[, func := factor(func, levels = c("truncating/splice","missense","in-frame indel"))]
u[, role := factor(role, levels = c("LoF","Act","ambiguous"))]
u[, is_hcc := as.logical(is_hcc)]

## ---- Panel A: per-gene mutation map (unique variants), grouped by role ----
# order genes: LoF first (by max impact then AF), then Act, then ambiguous
role_rank <- c(LoF = 1, Act = 2, ambiguous = 3)
u[, rrank := role_rank[as.character(role)]]
gene_ord <- u[, .(rrank = min(rrank), hi = any(impact == "HIGH"), maf = max(AF, na.rm = TRUE)),
              by = human_symbol][order(rrank, -hi, -maf)]
u[, gene := factor(human_symbol, levels = gene_ord$human_symbol)]

pA <- ggplot(u, aes(x = AF, y = gene)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "grey60") +
  geom_segment(aes(x = 0, xend = AF, y = gene, yend = gene, color = role), linewidth = 0.4) +
  geom_point(aes(color = role, shape = func, size = is_hcc), stroke = 0.4) +
  scale_color_manual(values = c(LoF = HEX$blue, Act = HEX$verm, ambiguous = HEX$grey),
                     labels = c(LoF = "Tumour suppressor (LoF)", Act = "Oncogene (GoF)", ambiguous = "Ambiguous")) +
  scale_shape_manual(values = c("truncating/splice" = 17, "missense" = 16, "in-frame indel" = 15)) +
  scale_size_manual(values = c(`TRUE` = 2.6, `FALSE` = 1.8), guide = "none") +
  scale_x_continuous(limits = c(0,1), expand = expansion(mult = c(0.02,0.03))) +
  labs(x = "Variant allele fraction", y = NULL, color = "Gene role",
       shape = "Consequence", title = "A") +
  theme_pub + theme(axis.text.y = element_text(size = 6.5, face = "italic"),
                    legend.position = "right")

## ---- Panel B: summary — number of mutated cancer genes by role x consequence ----
summ <- u[, .N, by = .(role, func)]
summ[, role := factor(role, levels = c("LoF","Act","ambiguous"))]
pB <- ggplot(summ, aes(x = role, y = N, fill = func)) +
  geom_col(position = "stack", width = 0.65, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("truncating/splice" = HEX$verm, "missense" = HEX$sky,
                               "in-frame indel" = HEX$green)) +
  scale_x_discrete(labels = c(LoF = "Tumour\nsuppressor", Act = "Oncogene", ambiguous = "Ambiguous")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Unique coding variants", fill = "Consequence", title = "B") +
  theme_pub + theme(legend.position = "right")

fig <- pA + pB + plot_layout(widths = c(1.7, 1))
svglite(file.path(OUT, "fig6_lof_gof_mutations.svg"), width = 10, height = 7, fix_text_size = FALSE)
print(fig); dev.off()
cat("wrote fig6_lof_gof_mutations.svg\n")
cat("unique variants by role:\n"); print(u[, .N, by = role])
cat("truncating/splice by role:\n"); print(u[func=="truncating/splice", .N, by = role])
