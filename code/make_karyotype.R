#!/usr/bin/env Rscript
# Molecular karyotype of the HeDe rat HCC cell line, derived from DNAcopy CBS
# copy-number segmentation (20 kb bins). Produces:
#   - karyotype_by_chromosome.tsv / .csv : per-chromosome copy-number summary
#   - report_figures/karyotype_ideogram.png / .svg : state-coloured ideogram
# Run from /workspace/hede_followup.  Preferred font Liberation Sans; Okabe-Ito palette.
#
# ---------------------------------------------------------------------------
# PROVENANCE NOTE (restored file)
# The original of this script was lost when its sandbox was recreated. It has been
# restored from (a) verbatim fragments preserved in the session transcript -- the
# header above, the library call, the two input-reading lines, the summarise
# pipeline and the "absent" call branch -- and (b) the surviving outputs.
# It is therefore a RECONSTRUCTION VERIFIED BY OUTPUT, not a byte-recovered
# original: re-running it regenerates karyotype_by_chromosome.tsv byte-for-byte
# (md5 10a052164fa27895eefa5d6a27755050) and .csv (md5 2af60818254bab8df3755e07d57b2e7a),
# and an ideogram matching the published figure.
# One parameter is under-determined by the data: the partial-gain/partial-loss
# fraction cut-off. Any value in (0.074, 0.172] reproduces the published calls
# exactly, because no chromosome has an altered fraction in that interval;
# 0.15 is used below. Whole-chromosome calls (>= 0.50) are unambiguous.
# ---------------------------------------------------------------------------

suppressMessages({library(ggplot2); library(dplyr); library(readr); library(svglite)})

WD  <- "/workspace/hede_followup"
DEP <- "/mnt/results/deposition/zenodo/HeDe_PacBio_GRCr8_dataset"
FIG <- file.path(WD, "report_figures"); dir.create(FIG, showWarnings = FALSE)
FONT <- "Liberation Sans"

# Prefer the deposited bundle so this script also runs from the Zenodo archive as-is;
# fall back to the original session paths that were used to make the published figure.
first_existing <- function(...) { p <- c(...); p[file.exists(p)][1] }
SEG_PATH <- first_existing(file.path(DEP, "cnv/HeDe_cnv_segments.tsv"),
                           "/mnt/shared-workspace/shared/HeDe_cnv_segments.tsv")
LEN_PATH <- first_existing(file.path(DEP, "tables/GRCr8_chromosome_lengths.tsv"),
                           "/mnt/shared-workspace/shared/rn_GRCr8.fa.fai")

# NOTE: HeDe_cnv_segments.tsv is comma-separated despite the .tsv extension.
seg <- read_csv(SEG_PATH, show_col_types = FALSE)
cat("Segments read:", nrow(seg), "\nUnique states:\n"); print(table(seg$state))

chr_order <- c(paste0("chr", 1:20), "chrX", "chrY")

# Chromosome lengths: GRCr8 (GCF_036323735.1). NC_0860XX.1 -> chr(XX-18); 39 -> X, 40 -> Y.
if (grepl("GRCr8_chromosome_lengths", LEN_PATH)) {
  len <- read_tsv(LEN_PATH, show_col_types = FALSE) %>%
    select(chr, len = length)
} else {
  fai <- read_tsv(LEN_PATH, col_names = FALSE, show_col_types = FALSE)
  n   <- suppressWarnings(as.integer(sub("^NC_0860([0-9]+)\\.1$", "\\1", fai$X1)))
  len <- tibble(
    chr = ifelse(!is.na(n) & n >= 19 & n <= 38, paste0("chr", n - 18),
          ifelse(!is.na(n) & n == 39, "chrX",
          ifelse(!is.na(n) & n == 40, "chrY", NA_character_))),
    len = fai$X2) %>% filter(!is.na(chr))
}
len <- len %>% filter(chr %in% chr_order)
cat("Chromosomes with lengths:", nrow(len), "\n")

seg <- seg %>% filter(chr %in% chr_order)

summ <- seg %>%
  group_by(chr) %>%
  summarise(
    n_seg      = n(),
    width_gain = sum(width[state == "gain"]),
    width_loss = sum(width[state == "loss"]),
    width_neut = sum(width[state == "neutral"]),
    wlog2r     = round(weighted.mean(log2r, width), 3),
    .groups    = "drop") %>%
  right_join(len, by = "chr") %>%
  mutate(
    n_seg      = ifelse(is.na(n_seg), 0, n_seg),
    width_gain = ifelse(is.na(width_gain), 0, width_gain),
    width_loss = ifelse(is.na(width_loss), 0, width_loss),
    width_neut = ifelse(is.na(width_neut), 0, width_neut),
    covered    = width_gain + width_loss + width_neut,
    frac_gain  = round(width_gain / len, 3),
    frac_loss  = round(width_loss / len, 3),
    frac_neut  = round(width_neut / len, 3),
    call = case_when(
      n_seg == 0        ~ "absent",
      frac_gain >= 0.50 ~ "gain",
      frac_loss >= 0.50 ~ "loss",
      frac_gain >= 0.15 ~ "partial gain",
      frac_loss >= 0.15 ~ "partial loss",
      TRUE              ~ "neutral")) %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  arrange(chr) %>%
  mutate(chr = as.character(chr)) %>%
  select(chr, n_seg, width_gain, width_loss, width_neut, wlog2r,
         len, covered, frac_gain, frac_loss, frac_neut, call)

cat("\n=== Per-chromosome molecular karyotype ===\n")
print(as.data.frame(summ))

write_tsv(summ, file.path(WD, "karyotype_by_chromosome.tsv"))
write.csv(summ, file.path(WD, "karyotype_by_chromosome.csv"), row.names = FALSE)

# genome-wide state totals (reconcile with earlier 1,078-segment summary)
cat("\n=== Genome-wide state totals (all segments) ===\n")
print(seg %>% count(state))

## ---------------------------------------------------------------- ideogram
OI_gain <- "#D55E00"   # Okabe-Ito vermillion
OI_loss <- "#0072B2"   # Okabe-Ito blue
OI_neut <- "#D9D9D9"   # neutral grey
# Bar half-height fitted to the surviving published PNG (bar 65 px, pitch 75 px at
# 8.5 x 6.5 in / 300 dpi); with zero y-expansion this reproduces its geometry exactly.
HH <- 0.42

pdat <- seg %>%
  mutate(chr = factor(chr, levels = rev(chr_order)),
         state = factor(state, levels = c("gain", "neutral", "loss")))
back <- len %>% mutate(chr = factor(chr, levels = rev(chr_order)))

p <- ggplot() +
  geom_rect(data = back,
            aes(xmin = 0, xmax = len / 1e6,
                ymin = as.numeric(chr) - HH, ymax = as.numeric(chr) + HH),
            fill = "#F4F4F4", colour = "#D9D9D9", linewidth = 0.25) +
  geom_rect(data = pdat,
            aes(xmin = start / 1e6, xmax = end / 1e6,
                ymin = as.numeric(chr) - HH, ymax = as.numeric(chr) + HH,
                fill = state),
            colour = NA) +
  scale_fill_manual(values = c(gain = OI_gain, neutral = OI_neut, loss = OI_loss),
                    name = "Copy-number\nstate") +
  scale_y_continuous(breaks = seq_along(chr_order),
                     labels = rev(chr_order), expand = expansion(add = 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  labs(title = "HeDe molecular karyotype",
       subtitle = "DNAcopy CBS copy-number segments (20 kb bins) per GRCr8 chromosome",
       x = "Position (Mb)", y = NULL) +
  theme_minimal(base_size = 12, base_family = FONT) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "#333333", linewidth = 0.5),
        axis.text = element_text(colour = "#4D4D4D"),
        axis.title = element_text(colour = "#4D4D4D"),
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 12, colour = "#4D4D4D"),
        legend.title = element_text(colour = "#4D4D4D"),
        legend.text = element_text(colour = "#4D4D4D"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

svglite(file.path(FIG, "karyotype_ideogram.svg"), width = 8.5, height = 6.5)
print(p); invisible(dev.off())
png(file.path(FIG, "karyotype_ideogram.png"), width = 8.5, height = 6.5,
    units = "in", res = 300)
print(p); invisible(dev.off())

cat("\nWrote:", file.path(WD, "karyotype_by_chromosome.tsv"), "\n")
cat("Wrote:", file.path(FIG, "karyotype_ideogram.png"), "/ .svg\n")
