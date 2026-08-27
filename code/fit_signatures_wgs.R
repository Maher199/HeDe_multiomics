#!/usr/bin/env Rscript
# Phase 7: refit COSMIC v3.3 SBS96 signatures on the WGS somatic spectrum and
# compare with the PacBio tumor-only fit.
#
# Method validated to reproduce the deposited PacBio fit EXACTLY (all 14
# signatures, max abs diff 0 of 2,231,053): NNLS (nnls package) on proportions
# against the cosmicsig COSMIC_v3.3$signature$GRCh37$SBS96 matrix, with the
# cosmicsig compact channel names mapped as [5'][ref][3'][alt] -> 5'[ref>alt]3'.
#
# Usage: Rscript fit_signatures_wgs.R <wgs.96spectrum.csv> <outdir>
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages({library(nnls); library(cosmicsig); library(readr); library(dplyr); library(tidyr)})

args <- commandArgs(trailingOnly = TRUE)
spec_f <- args[1]; outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

PAC_TGT <- "/mnt/results/deposition/zenodo/HeDe_PacBio_GRCr8_dataset/tables/signature_contributions_targeted.csv"

sigs <- COSMIC_v3.3$signature$GRCh37$SBS96
compact_to_long <- function(x) paste0(substr(x, 1, 1), "[", substr(x, 2, 2), ">",
                                      substr(x, 4, 4), "]", substr(x, 3, 3))

spec <- read_csv(spec_f, show_col_types = FALSE)
stopifnot(nrow(spec) == 96)
v <- spec$count[match(compact_to_long(rownames(sigs)), spec$channel)]
if (any(is.na(v))) stop("channel-name mismatch between spectrum and cosmicsig matrix")
cat("total somatic SNVs in spectrum:", sum(v), "\n")

panel14 <- c("SBS1","SBS4","SBS5","SBS8","SBS11","SBS19","SBS22","SBS24",
             "SBS30","SBS37","SBS39","SBS40","SBS42","SBS44")
stopifnot(all(panel14 %in% colnames(sigs)))

fit_one <- function(sig_mat) {
  fit <- nnls(as.matrix(sig_mat), v / sum(v))
  n <- round(fit$x * sum(v))
  data.frame(signature = colnames(sig_mat), n = n,
             frac = round(n / sum(v), 4), row.names = NULL) %>%
    arrange(desc(frac))
}

tgt <- fit_one(sigs[, panel14])
unr <- fit_one(sigs)

write_csv(tgt, file.path(outdir, "HeDe_WGS_signature_targeted.csv"))
write_csv(unr, file.path(outdir, "HeDe_WGS_signature_contributions.csv"))

pac_t <- read_csv(PAC_TGT, show_col_types = FALSE)
cmp <- full_join(pac_t %>% rename(n_pacbio = n, frac_pacbio = frac),
                 tgt %>% rename(n_wgs = n, frac_wgs = frac),
                 by = "signature") %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  arrange(desc(frac_wgs))
write_csv(cmp, file.path(outdir, "signature_fit_pacbio_vs_wgs.csv"))

cat("\n=== targeted fit: WGS somatic vs PacBio tumor-only ===\n")
print(as.data.frame(cmp), row.names = FALSE)

# cosine similarity between the two observed spectra (PacBio novel vs WGS somatic)
pac_spec <- read_csv("/mnt/results/deposition/zenodo/HeDe_PacBio_GRCr8_dataset/tables/mutational_spectrum_96.csv",
                     show_col_types = FALSE)
v_pac <- pac_spec$count[match(spec$channel, pac_spec$channel)]
cos_sim <- sum(v * v_pac) / sqrt(sum(v^2) * sum(v_pac^2))
cat(sprintf("\ncosine similarity, WGS-somatic vs PacBio-novel 96-spectra: %.4f\n", cos_sim))
cat("DONE\n")
