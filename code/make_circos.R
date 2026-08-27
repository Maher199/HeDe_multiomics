#!/usr/bin/env Rscript
# Genome-wide Circos-style integration plot: SVs + CNV + coding mutations + CGI hypermethylation.
suppressMessages({library(data.table); library(circlize); library(svglite)})

WD  <- "/workspace/hede_followup"
OUT <- "/workspace/hede_followup/figures"; dir.create(OUT, showWarnings = FALSE)
FONT <- "Liberation Sans"
HEX <- list(grey="#999999", orange="#E69F00", sky="#56B4E9", green="#009E73",
            yellow="#F0E442", blue="#0072B2", verm="#D55E00", pink="#CC79A7",
            black="#000000", lgrey="#E0E0E0")

## ---- chromosomes (chr1-20, X, Y) ----
fai <- fread("/mnt/shared-workspace/shared/rn_GRCr8.fa.fai", col.names = c("chrom","len","a","b","c"))
fai <- fai[grepl("^NC_0860(19|2[0-9]|3[0-9]|40)", chrom)]
ncnum <- as.integer(sub("NC_0860", "", fai$chrom))
fai[, chr := ifelse(ncnum <= 38, as.character(ncnum - 18), ifelse(ncnum == 39, "X", "Y"))]
fai[, chr := factor(chr, levels = c(1:20, "X", "Y"))]
setorder(fai, chr)
chroms <- as.character(fai$chrom)
chr_label <- setNames(as.character(fai$chr), fai$chrom)
genome_len <- setNames(fai$len, fai$chrom)

## ---- CNV (chrN -> NC) ----
cnv <- fread("/mnt/shared-workspace/shared/HeDe_cnv_segments.tsv")
chr2nc <- setNames(fai$chrom, as.character(fai$chr))
cnv[, chrom := chr2nc[sub("^chr", "", chr)]]
cnv <- cnv[!is.na(chrom)]
cnv[, log2r_c := pmax(-3, pmin(3, as.numeric(log2r)))]

## ---- SVs ----
sv <- fread(file.path(WD, "variants/sv_parsed.tsv"), header = FALSE,
            col.names = c("chrom","pos","svtype","svlen","chr2","end2","support"))
sv <- sv[chrom %in% chroms]
sv[, svlen := suppressWarnings(as.numeric(svlen))]
sv[, support := suppressWarnings(as.numeric(support))]
BIN <- 1e6
sv[, bin := floor(pos / BIN) * BIN]
svdens <- sv[, .N, by = .(chrom, bin)]
bigsv <- sv[svtype %in% c("DEL","DUP","INV") & !is.na(svlen) & abs(svlen) >= 1e6]

## ---- BND translocation mates (interchromosomal, SUPPORT>=10) ----
bnd <- fread(file.path(WD, "variants/bnd_mates.tsv"), header = FALSE,
             col.names = c("chrom","pos","chr2","mate_chrom","mate_pos","support"))
bnd <- bnd[!is.na(chrom) & !is.na(mate_chrom) & !is.na(pos) & !is.na(mate_pos)]
bnd <- bnd[chrom %in% chroms & mate_chrom %in% chroms & chrom != mate_chrom & support >= 30]
bnd[, k1 := pmin(paste(chrom, pos), paste(mate_chrom, mate_pos))]
bnd[, k2 := pmax(paste(chrom, pos), paste(mate_chrom, mate_pos))]
bnd <- unique(bnd, by = c("k1","k2"))

## ---- coding mutations ----
cm  <- fread(file.path(WD, "variants/coding_mutations_all.tsv"))[CHROM %in% chroms]
cmu <- unique(cm, by = c("CHROM","POS","REF","ALT"))
cc  <- fread(file.path(WD, "variants/cancer_coding_mutations.tsv"))
ccu <- unique(cc, by = c("human_symbol","CHROM","POS"))[CHROM %in% chroms]

## ---- CGI hypermethylation ----
cgi <- fread(file.path(WD, "methyl/cgi_5mC.tsv"))
cgih <- cgi[meth_state == "hypermethylated" & chrom %in% chroms]

## ---- driver loci (coding + epigenetic-silencing hits; gene midpoints) ----
drivers <- data.table(
  label = c("Cdkn2a/b", "Tp53", "SMARCA4", "LRP1B", "RPS6KA3"),
  chrom = c("NC_086023.1", "NC_086028.1", "NC_086026.1", "NC_086021.1", "NC_086039.1"),
  pos   = c(109100000, 54807961, 28486720, 46064574, 39378971))
drivers <- drivers[chrom %in% chroms]

cat("CNV", nrow(cnv), "| SVdens", nrow(svdens), "| largeSV", nrow(bigsv),
    "| BNDarcs", nrow(bnd), "| coding", nrow(cmu), "| cancer", nrow(ccu),
    "| CGIhyper", nrow(cgih), "| drivers", nrow(drivers), "\n")

## ================= render =================
svgfile <- file.path(OUT, "fig5_genomewide_circos.svg")
svglite(svgfile, width = 9.5, height = 9.5, fix_text_size = FALSE)
par(family = FONT, mar = c(2.5,2.5,2.5,2.5), xpd = NA)

circos.clear()
circos.par(start.degree = 90, gap.degree = 1.5, cell.padding = c(0,0,0,0),
           track.margin = c(0.008, 0.008), points.overflow.warning = FALSE,
           canvas.xlim = c(-1.15, 1.15), canvas.ylim = c(-1.15, 1.15))
xlim_mat <- matrix(c(rep(0, length(chroms)), as.numeric(genome_len[chroms])), ncol = 2)
circos.initialize(factors = chroms, xlim = xlim_mat)

## Track 1: ideogram + labels
circos.track(ylim = c(0,1), track.height = 0.07, bg.border = NA,
             panel.fun = function(x, y) {
               chr <- CELL_META$sector.index
               circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                           col = "#F4F4F4", border = "grey60", lwd = 0.3)
               circos.text(CELL_META$xcenter, 0.5, chr_label[[chr]],
                           cex = 0.6, facing = "inside", niceFacing = TRUE, font = 2)
             })

## Track 2: CNV log2r
circos.track(ylim = c(-3, 3), track.height = 0.13, bg.border = NA,
             panel.fun = function(x, y) {
               chr <- CELL_META$sector.index
               d <- cnv[chrom == chr]; if (!nrow(d)) return()
               circos.lines(CELL_META$cell.xlim, c(0,0), col = "grey85", lwd = 0.3)
               col <- ifelse(d$log2r_c > 0.25, HEX$verm, ifelse(d$log2r_c < -0.25, HEX$blue, HEX$lgrey))
               circos.points((d$start+d$end)/2, d$log2r_c, col = col, pch = 16, cex = 0.3)
             })

## Track 3: SV density histogram
ymax <- max(svdens$N)
circos.track(ylim = c(0, ymax), track.height = 0.10, bg.border = NA,
             panel.fun = function(x, y) {
               chr <- CELL_META$sector.index
               d <- svdens[chrom == chr]; if (!nrow(d)) return()
               circos.barplot(d$N, d$bin, col = HEX$sky, border = NA, bar_width = BIN)
             })

## Track 4: large SVs (>=1 Mb)
circos.track(ylim = c(0,1), track.height = 0.05, bg.border = NA,
             panel.fun = function(x, y) {
               chr <- CELL_META$sector.index
               d <- bigsv[chrom == chr]; if (!nrow(d)) return()
               col <- ifelse(d$svtype=="DEL", HEX$blue, ifelse(d$svtype=="DUP", HEX$verm, HEX$green))
               circos.rect(d$pos, 0, d$pos + abs(d$svlen), 1, col = col, border = NA)
             })

## Track 5: coding mutations
circos.track(ylim = c(0,1), track.height = 0.08, bg.border = NA,
             panel.fun = function(x, y) {
               chr <- CELL_META$sector.index
               d <- cmu[CHROM == chr]
               if (nrow(d)) circos.points(d$POS, rep(0.35, nrow(d)), col = "#00000022", pch = 16, cex = 0.15)
               dc <- ccu[CHROM == chr]
               if (nrow(dc)) {
                 col <- ifelse(dc$impact == "HIGH", HEX$verm, HEX$orange)
                 circos.points(dc$POS, rep(0.75, nrow(dc)), col = col, pch = 17, cex = 0.55)
               }
             })

## Track 6: CGI hypermethylation
circos.track(ylim = c(0,1), track.height = 0.05, bg.border = NA,
             panel.fun = function(x, y) {
               chr <- CELL_META$sector.index
               d <- cgih[chrom == chr]; if (!nrow(d)) return()
               circos.points((d$start+d$end)/2, rep(0.5, nrow(d)), col = HEX$green, pch = 16, cex = 0.2)
             })

## Center: translocation arcs
if (nrow(bnd)) {
  for (i in seq_len(nrow(bnd))) {
    r <- bnd[i]
    c1 <- as.character(r$chrom); c2 <- as.character(r$mate_chrom)
    tryCatch(circos.link(c1, r$pos, c2, r$mate_pos, col = "#0072B2A0", border = NA, lwd = 0.8, rou = 0.5),
             error = function(e) NULL)
  }
}

## Driver-locus markers + labels (ideogram, pointing outward, clear of the ring)
for (i in seq_len(nrow(drivers))) {
  d <- drivers[i]
  circos.points(d$pos, 0.5, sector.index = d$chrom, track.index = 1,
                pch = 25, col = HEX$black, bg = HEX$yellow, cex = 0.85)
  # place label outside the ideogram, radially oriented and auto-flipped for readability
  circos.text(d$pos, 2.4, d$label, sector.index = d$chrom, track.index = 1,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
              cex = 0.7, col = HEX$black, font = 3)
}

circos.clear()

## Legend (base graphics, placed in lower-left inside the canvas)
par(xpd = NA)
legend(x = -1.13, y = 0.55, bty = "n", cex = 0.66, ncol = 1, xjust = 0, yjust = 0,
       legend = c("CNV gain", "CNV loss", "SV density", "Large SV (DEL/DUP/INV)",
                  "Coding mutation", "Cancer gene (HIGH/MOD)", "CGI hypermethylated", "Translocation"),
       pch = c(16,16,15,15,16,17,16,NA), lty = c(NA,NA,NA,NA,NA,NA,NA,1),
       col = c(HEX$verm, HEX$blue, HEX$sky, HEX$green, "grey40", HEX$orange, HEX$green, HEX$blue),
       pt.cex = c(1,1,1.2,1.2,0.8,1,0.8,NA), lwd = c(NA,NA,NA,NA,NA,NA,NA,1.5))

dev.off()
cat("wrote", svgfile, "\n")
