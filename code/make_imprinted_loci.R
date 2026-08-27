#!/usr/bin/env Rscript
# Define rat imprinted-loci regions (germline DMRs / ICRs) for methylation profiling.
# Strategy: promoter DMR (TSS +/- 2 kb, strand-aware) for each imprinted gene, plus the two
# intergenic ICRs (H19 ICR; Dlk1-Gtl2 IG-DMR) defined from conserved synteny. Coordinates are
# approximate for the intergenic ICRs (documented as such); promoter DMRs are exact from the GFF.
suppressMessages({library(data.table)})
genes <- fread("/workspace/hede_followup/ref/GRCr8_genes.bed",
               col.names = c("chrom","start","end","strand","symbol"))
imp <- c("Igf2","H19","Dlk1","Meg3","Rtl1","Peg3","Snrpn","Kcnq1","Grb10","Igf2r",
         "Mest","Peg10","Nnat","Impact","Plagl1","Cdkn1c","Phlda2","Gnas","Rasgrf1",
         "Dcn","Slc22a18","Ascl2","Tspan32","Cd81","Trpm5","Trappc9","Ppp1r9a")
g <- genes[symbol %in% imp]
HALF <- 2000L
# strand-aware TSS: + strand -> start; - strand -> end
g[, tss := ifelse(strand == "+", start, end)]
g[, rstart := pmax(0L, tss - HALF)]
g[, rend   := tss + HALF]
prom <- g[, .(chrom, start = rstart, end = rend, region = paste0(symbol, "_promDMR"), class = "promoter_DMR")]

# Intergenic ICRs (approximate, conserved synteny with mouse/human)
icr <- data.table(
  chrom = c("NC_086019.1", "NC_086024.1"),
  start = c(207163000, 134256000),   # H19 ICR (H19-Igf2 intergenic); IG-DMR (Dlk1-Meg3 intergenic)
  end   = c(207168000, 134261000),
  region = c("H19_ICR", "Dlk1_Meg3_IGDMR"),
  class  = c("intergenic_ICR", "intergenic_ICR")
)
out <- rbind(prom, icr)
setorder(out, chrom, start)
fwrite(out, "/workspace/hede_followup/ref/imprinted_loci.bed", sep = "\t",
       col.names = FALSE)
cat("imprinted regions written:", nrow(out), "\n")
print(out)
