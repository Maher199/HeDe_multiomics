#!/usr/bin/env Rscript
# Phase 6 parser: classify SURVIVOR-merged SVs by support vector (01 = Sniffles-only,
# 10 = Manta-only, 11 = both) and summarise concordance by SV type and size.
# Usage: Rscript sv_concordance.R <merged_survivor.vcf> <outdir>
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages({library(readr); library(dplyr); library(stringr); library(tidyr)})

args <- commandArgs(trailingOnly = TRUE)
vcf <- args[1]; outdir <- args[2]

lines <- readLines(vcf)
hdr <- grep("^##", lines)
dat <- read_tsv(I(lines[-hdr]), col_names = FALSE, show_col_types = FALSE,
                col_types = cols(.default = "c"))
colnames(dat)[1:9] <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

sv <- dat %>%
  mutate(
    SVTYPE = str_extract(INFO, "SVTYPE=[A-Z]+") %>% str_remove("SVTYPE="),
    SVLEN  = as.numeric(str_extract(INFO, "SVLEN=-?[0-9]+") %>% str_remove("SVLEN=")),
    SUPP   = str_extract(INFO, "SUPP=[0-9]+") %>% str_remove("SUPP="),
    SUPP_VEC = str_extract(INFO, "SUPP_VEC=[01]+") %>% str_remove("SUPP_VEC="),
    support = case_when(SUPP_VEC == "11" ~ "both",
                        SUPP_VEC == "10" ~ "manta_only",
                        SUPP_VEC == "01" ~ "sniffles_only",
                        TRUE ~ "other"),
    size_class = case_when(abs(SVLEN) < 1000 ~ "50-1k",
                           abs(SVLEN) < 10000 ~ "1k-10k",
                           abs(SVLEN) < 100000 ~ "10k-100k",
                           TRUE ~ ">=100k"))

write_tsv(sv %>% select(CHROM, POS, SVTYPE, SVLEN, support, SUPP_VEC),
          file.path(outdir, "sv_concordance_merged.tsv"))

by_type <- sv %>% count(SVTYPE, support) %>% pivot_wider(names_from = support, values_from = n, values_fill = 0)
by_size <- sv %>% count(size_class, support) %>% pivot_wider(names_from = support, values_from = n, values_fill = 0)
write_tsv(by_type, file.path(outdir, "sv_concordance_by_type.tsv"))
write_tsv(by_size, file.path(outdir, "sv_concordance_by_size.tsv"))

tot <- sv %>% count(support)
cat("=== SV concordance (SURVIVOR, 500 bp, type-matched) ===\n")
print(as.data.frame(tot), row.names = FALSE)
cat("\nby type:\n");  print(as.data.frame(by_type), row.names = FALSE)
cat("\nby size:\n");  print(as.data.frame(by_size), row.names = FALSE)

m_total <- sum(sv$support %in% c("manta_only", "both"))
s_total <- sum(sv$support %in% c("sniffles_only", "both"))
both    <- sum(sv$support == "both")
cat(sprintf("\nManta SVs: %d | Sniffles SVs: %d | shared: %d\n", m_total, s_total, both))
cat(sprintf("concordance: %.1f%% of Manta, %.1f%% of Sniffles\n",
            100 * both / m_total, 100 * both / s_total))
