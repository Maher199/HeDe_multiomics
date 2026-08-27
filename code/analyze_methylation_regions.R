#!/usr/bin/env Rscript
# Region-level 5mC profiling for HeDe (single-sample aberrant-methylation framework).
# Aggregates per-CpG 5mC over CpG islands, retrotransposons (LINE/SINE/LTR), and imprinted
# loci, using coverage-weighted means. Chunked by chromosome to bound memory.
suppressMessages({library(data.table)})
setDTthreads(1)
ref <- "/workspace/hede_followup/ref"
METH <- "/mnt/shared-workspace/shared/methyl/HeDe_cpg_5mC_filtered.tsv"
OUTD <- "/workspace/hede_followup/methyl"; dir.create(OUTD, showWarnings = FALSE)
MINCOV <- 10L

cgi <- fread(file.path(ref, "GRCr8_cgi.bed"), col.names = c("chrom","start","end","ncpg","gcpct","oe"))
cgi[, region_id := paste0("CGI_", .I)]; cgi[, class := "CpG_island"]

# RepeatMasker .out has 15 fields for '+' strand and 16 for '-' strand (extra paren col).
# Read without names, fill to 16, then select positionally:
# V5=chrom V6=begin V7=end V9=strand V10=repeat V11=class/family (constant across strand)
rm0 <- fread(file.path(ref, "GCF_036323735.1_GRCr8_rm.out.gz"), skip = 3, fill = TRUE,
             header = FALSE, sep = " ", quote = "")
setnames(rm0, paste0("V", seq_len(ncol(rm0))))
rm <- rm0[grepl("^(LINE|SINE|LTR)/", V11)]
rm[, `:=`(chrom = V5, begin = as.integer(V6), end = as.integer(V7),
          classfam = V11, repname = V10)]
rm[, rclass := sub("/.*", "", classfam)]
rm[, start := begin - 1L]
repeats <- rm[, .(chrom, start, end, rclass, classfam, repname)]
repeats <- repeats[!is.na(start) & !is.na(end) & end > start]
repeats[, region_id := paste0("REP_", .I)]
cat("retrotransposon repeats (LINE/SINE/LTR):", nrow(repeats), "\n")

imp <- fread(file.path(ref, "imprinted_loci.bed"),
             col.names = c("chrom","start","end","region_id","class"))

cat("loading methylation ...\n")
m <- fread(METH, select = c("chrom","start","end","Nvalid","frac"))
m <- m[Nvalid >= MINCOV]
m[, `:=`(start = as.integer(start), end = as.integer(end), Nvalid = as.integer(Nvalid))]
setkey(m, chrom, start, end)
cat("CpGs with Nvalid>=", MINCOV, ":", nrow(m), "\n")
genome_mean <- m[, weighted.mean(frac, Nvalid)]
cat("genome-wide coverage-weighted mean 5mC:", round(genome_mean, 4), "\n")

aggregate_regions <- function(regions) {
  out <- vector("list", length(unique(m$chrom))); i <- 0
  for (ch in unique(m$chrom)) {
    r <- regions[chrom == ch]; if (nrow(r) == 0) next
    mm <- m[chrom == ch]; setkey(r, start, end)
    ov <- foverlaps(mm, r, by.x = c("start","end"), by.y = c("start","end"), type = "any", nomatch = NULL)
    if (nrow(ov) == 0) next
    i <- i + 1
    out[[i]] <- ov[, .(n_cpg = .N, mean5mC = weighted.mean(frac, Nvalid),
                       total_cov = sum(Nvalid)), by = region_id]
  }
  rbindlist(out)
}

cat("aggregating CGI ...\n")
cgi_stats <- aggregate_regions(cgi[, .(chrom, start, end, region_id)])
cgi_stats <- merge(cgi_stats, cgi[, .(region_id, chrom, start, end)], by = "region_id")
cgi_stats[, meth_state := fifelse(mean5mC >= 0.6, "hypermethylated",
                          fifelse(mean5mC <= 0.2, "hypomethylated", "intermediate"))]
fwrite(cgi_stats, file.path(OUTD, "cgi_5mC.tsv"), sep = "\t")
print(cgi_stats[, .N, by = meth_state])

cat("aggregating retrotransposons by class ...\n")
rep_list <- vector("list", length(unique(m$chrom))); i <- 0
for (ch in unique(m$chrom)) {
  r <- repeats[chrom == ch]; if (nrow(r) == 0) next
  mm <- m[chrom == ch]; setkey(r, start, end)
  ov <- foverlaps(mm, r, by.x = c("start","end"), by.y = c("start","end"), type = "any", nomatch = NULL)
  if (nrow(ov) == 0) next
  i <- i + 1
  rep_list[[i]] <- ov[, .(n_cpg = .N, total_cov = sum(Nvalid),
                          wsum = sum(frac * Nvalid)), by = .(rclass, classfam)]
}
rep_class <- rbindlist(rep_list)[, .(n_cpg = sum(n_cpg), total_cov = sum(total_cov),
                                     wsum = sum(wsum)), by = .(rclass, classfam)]
rep_class[, mean5mC := wsum / total_cov][, wsum := NULL]
rep_class[, genome_mean := genome_mean][, delta_vs_genome := mean5mC - genome_mean]
setorder(rep_class, rclass, -n_cpg)
fwrite(rep_class, file.path(OUTD, "retrotransposon_5mC.tsv"), sep = "\t")
print(rep_class)

cat("aggregating imprinted loci ...\n")
imp_stats <- aggregate_regions(imp[, .(chrom, start, end, region_id)])
imp_stats <- merge(imp_stats, imp[, .(region_id, chrom, start, end, class)], by = "region_id")
imp_stats[, loi_state := fifelse(mean5mC < 0.35, "loss(hypo)",
                          fifelse(mean5mC > 0.65, "gain(hyper)", "intact(~0.5)"))]
fwrite(imp_stats, file.path(OUTD, "imprinted_5mC.tsv"), sep = "\t")
print(imp_stats[, .(region_id, class, n_cpg, mean5mC, loi_state)])
cat("DONE\n")
