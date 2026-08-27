#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Recover the exact per-read length distribution of a PacBio HiFi BAM from its
# companion .pbi index, and compute yield / N50 / N90 / percentiles.
#
# This script exists because the read-length statistics of this dataset must be
# auditable from primary data. An earlier version of Figure 1 and of Table 2 in
# the Data Descriptor reported a read N50 of 23,526 bp and a mean of 23,243 bp.
# Those values came from the first 500,000 reads in BAM order, which is a
# length-biased head-of-file subsample (early ZMWs skew long) rather than a
# random sample. Parsing the complete index gives N50 = 23,294 bp and
# mean = 23,020 bp over all 5,565,830 reads; those are the published values.
#
# Usage:
#   Rscript read_lengths_from_pbi.R <in.hifi_reads.bam.pbi> <outdir>
#
# Outputs (written to <outdir>):
#   read_lengths_full.txt.gz  one integer read length per line, BAM order
#   read_length_summary.tsv   yield, mean, median, N50, N90, min/max, percentiles
#
# The .pbi is a gzip-compressed binary file. Layout (pbbam specification):
#   header      magic "PBI\1" (4 B) | version uint32 (4 B) | pbiFlags uint16 (2 B)
#               | nReads uint32 (4 B) | reserved (18 B)                  = 32 B
#   BasicData   columnar arrays, each of length nReads, in this order:
#               rgId int32, qStart int32, qEnd int32, holeNumber int32,
#               readQual float32, ctxtFlag uint8, fileOffset int64  = 29 B/read
#   BarcodeData bcForward int16, bcReverse int16, bcQual int8        =  5 B/read
#               (present only when pbiFlags has the BARCODE bit set)
# Read length = qEnd - qStart. Only rgId, qStart and qEnd need to be read.
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript read_lengths_from_pbi.R <in.bam.pbi> <outdir>")
}
pbi    <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

con <- gzfile(pbi, "rb")

magic <- readBin(con, "raw", n = 4)
if (!identical(as.integer(magic), c(0x50L, 0x42L, 0x49L, 0x01L))) {
  close(con); stop("not a PBI file: magic bytes are not \"PBI\\1\"")
}
version  <- readBin(con, "integer", n = 1, size = 4, endian = "little")
pbiFlags <- readBin(con, "integer", n = 1, size = 2, signed = FALSE,
                    endian = "little")
nReads   <- readBin(con, "integer", n = 1, size = 4, endian = "little")
invisible(readBin(con, "raw", n = 18))   # reserved

has_bc <- bitwAnd(pbiFlags, 4L) != 0L    # BARCODE bit
cat(sprintf("PBI version=%d  pbiFlags=%d (barcode=%s)  nReads=%d\n",
            version, pbiFlags, has_bc, nReads))

# Self-check: the uncompressed size implied by the header must match the file.
# 32 B header + 29 B/read BasicData (+ 5 B/read BarcodeData when present).
bytes_per_read <- 29 + if (has_bc) 5 else 0
cat(sprintf("implied uncompressed size: %.0f B (%d B/read)\n",
            32 + bytes_per_read * as.numeric(nReads), bytes_per_read))

# BasicData is columnar: read rgId, then qStart, then qEnd, and stop there.
invisible(readBin(con, "integer", n = nReads, size = 4, endian = "little"))
qStart <- readBin(con, "integer", n = nReads, size = 4, endian = "little")
qEnd   <- readBin(con, "integer", n = nReads, size = 4, endian = "little")
close(con)

stopifnot(length(qStart) == nReads, length(qEnd) == nReads)

len <- as.numeric(qEnd) - as.numeric(qStart)
ok  <- is.finite(len) & len > 0
if (sum(ok) < nReads) {
  cat(sprintf("note: %d of %d records had non-positive length and were dropped\n",
              nReads - sum(ok), nReads))
}
len <- len[ok]

total <- sum(len)
l     <- sort(len, decreasing = TRUE)
cs    <- cumsum(l)
n50   <- l[which(cs >= total * 0.5)[1]]
n90   <- l[which(cs >= total * 0.9)[1]]
q     <- quantile(len, c(.01, .05, .25, .5, .75, .95, .99))

cat(sprintf("\n=== full-set read statistics (n = %d) ===\n", length(len)))
cat(sprintf("total yield   : %.0f bp (%.2f Gb)\n", total, total / 1e9))
cat(sprintf("mean length   : %.1f bp\n", mean(len)))
cat(sprintf("median length : %.0f bp\n", median(len)))
cat(sprintf("read N50      : %.0f bp (%.2f kb)\n", n50, n50 / 1000))
cat(sprintf("read N90      : %.0f bp (%.2f kb)\n", n90, n90 / 1000))
cat(sprintf("min / max     : %.0f / %.0f bp\n", min(len), max(len)))

gz <- gzfile(file.path(outdir, "read_lengths_full.txt.gz"), "wt")
writeLines(format(len, trim = TRUE, scientific = FALSE), gz)
close(gz)

summ <- data.frame(
  statistic = c("n_reads_indexed", "n_reads_usable", "total_yield_bp",
                "total_yield_Gb", "mean_length_bp", "median_length_bp",
                "read_N50_bp", "read_N90_bp", "min_length_bp", "max_length_bp",
                "p01_bp", "p05_bp", "p25_bp", "p75_bp", "p95_bp", "p99_bp"),
  value = c(nReads, length(len), total, round(total / 1e9, 3),
            round(mean(len), 1), median(len), n50, n90, min(len), max(len),
            round(as.numeric(q[1])), round(as.numeric(q[2])),
            round(as.numeric(q[3])), round(as.numeric(q[5])),
            round(as.numeric(q[6])), round(as.numeric(q[7]))),
  stringsAsFactors = FALSE)
write.table(summ, file.path(outdir, "read_length_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nwrote read_lengths_full.txt.gz (%d lines) and read_length_summary.tsv to %s\n",
            length(len), outdir))
