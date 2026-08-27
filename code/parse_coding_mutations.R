#!/usr/bin/env Rscript
# Analysis A: parse SnpEff-annotated novel PASS variants -> coding mutations -> cancer-gene overlap.
# Inputs : variants/novel_PASS.ann.vcf.gz, ref/cancer_genes_rat.tsv
# Outputs: variants/coding_mutations_all.tsv, variants/cancer_coding_mutations.tsv, summary to stdout
suppressMessages({ library(data.table) })

WD      <- "/workspace/hede_followup"
ann_vcf <- file.path(WD, "variants/novel_PASS.ann.vcf.gz")
cancer  <- fread(file.path(WD, "ref/cancer_genes_rat.tsv"))
cancer[, gene_key := tolower(rat_symbol)]

message("[1] Reading annotated VCF: ", ann_vcf)
vcf <- fread(cmd = paste("zcat", shQuote(ann_vcf), "| grep -v '^##'"),
             sep = "\t", header = TRUE, showProgress = FALSE)
setnames(vcf, 1, "CHROM")
setnames(vcf, names(vcf), sub("^#", "", names(vcf)))
message("    variants: ", nrow(vcf))

# --- parse FORMAT/SAMPLE for GT, DP, AD, AF (map by name; robust to field-count variation) ---
fmt <- vcf$FORMAT; smp <- vcf[[ncol(vcf)]]
fidx <- strsplit(fmt[1], ":", fixed = TRUE)[[1]]
sp   <- tstrsplit(smp, ":", fixed = TRUE, fill = NA_character_)
getf <- function(name) {
  i <- match(name, fidx)
  if (is.na(i) || length(sp) < i) return(rep(NA_character_, nrow(vcf)))
  sp[[i]]
}
vcf[, GT := getf("GT")]
vcf[, DP := suppressWarnings(as.integer(getf("DP")))]
vcf[, AF := suppressWarnings(as.numeric(sub(",.*$", "", getf("AF"))))]  # first (alt) AF
adr <- tstrsplit(getf("AD"), ",", fixed = TRUE, fill = NA_character_)
vcf[, AD_ref := suppressWarnings(as.integer(adr[[1]]))]
vcf[, AD_alt := suppressWarnings(as.integer(if (length(adr) >= 2) adr[[2]] else NA_character_))]
# fallback: derive AD from AF*DP where AD missing
bad <- is.na(vcf$AD_alt) & !is.na(vcf$AF) & !is.na(vcf$DP)
vcf$AD_alt[bad] <- round(vcf$AF[bad] * vcf$DP[bad])
vcf$AD_ref[bad] <- vcf$DP[bad] - vcf$AD_alt[bad]

# --- extract ANN subfield from INFO ---
info <- vcf$INFO
has  <- grepl("ANN=", info, fixed = TRUE)
message("    variants with ANN: ", sum(has))
ann  <- rep(NA_character_, length(info))
ann[has] <- sub(".*ANN=", "", info[has])          # strip everything up to ANN=
ann[has] <- sub(";.*$", "", ann[has])             # ANN value ends at ';' or EOL

# --- split ANN into individual annotations (long table) ---
if (sum(has) == 0) stop("No ANN annotations found in INFO field.")
ann_list <- strsplit(ann[has], ",", fixed = TRUE)
var_idx  <- rep(which(has), lengths(ann_list))
all_ann  <- unlist(ann_list, use.names = FALSE)
fp <- tstrsplit(all_ann, "|", fixed = TRUE, fill = "")
# ensure at least 11 fields WITHOUT truncating (pad only appends if a field is missing)
getf2 <- function(i) if (length(fp) >= i) fp[[i]] else rep("", length(all_ann))

long <- data.table(
  var_idx    = var_idx,
  allele     = getf2(1), annotation = getf2(2), impact = getf2(3),
  gene       = getf2(4), gene_id    = getf2(5), feature_type = getf2(6),
  feature_id = getf2(7), biotype    = getf2(8), rank = getf2(9),
  hgvs_c     = getf2(10), hgvs_p    = getf2(11)
)
long[, `:=`(CHROM = vcf$CHROM[var_idx], POS = vcf$POS[var_idx],
            REF = vcf$REF[var_idx], ALT = vcf$ALT[var_idx],
            GT = vcf$GT[var_idx], DP = vcf$DP[var_idx], AF = vcf$AF[var_idx],
            AD_ref = vcf$AD_ref[var_idx], AD_alt = vcf$AD_alt[var_idx])]
long[, gene_key := tolower(gene)]

# --- classify ---
CODING_IMP <- c("HIGH", "MODERATE")
long[, is_coding := impact %in% CODING_IMP]

# consequence group (collapse multi-term annotations like "missense_variant&splice_region")
consequence_group <- function(a) {
  a <- tolower(a)
  fifelse(grepl("frameshift", a), "frameshift",
  fifelse(grepl("stop_gained|stop_lost|start_lost", a), "nonsense/start-stop",
  fifelse(grepl("splice_acceptor|splice_donor|splice_region", a), "splice",
  fifelse(grepl("missense", a), "missense",
  fifelse(grepl("inframe_insertion|inframe_deletion|inframe", a), "inframe_indel",
  fifelse(grepl("synonymous", a), "synonymous",
  fifelse(grepl("utr", a), "UTR",
  fifelse(grepl("intron", a), "intron",
  fifelse(grepl("intergenic", a), "intergenic", "other")))))))))
}
long[, cons_group := consequence_group(annotation)]

# --- coding mutations table (one row per variant x coding-annotation) ---
coding <- long[is_coding == TRUE]
setorder(coding, CHROM, POS)
fwrite(coding[, .(CHROM, POS, REF, ALT, gene, gene_id, feature_id, biotype,
                  impact, annotation, cons_group, hgvs_c, hgvs_p,
                  GT, DP, AD_ref, AD_alt, AF)],
       file.path(WD, "variants/coding_mutations_all.tsv"), sep = "\t")
message("[2] coding (HIGH/MOD) annotations: ", nrow(coding),
        " | unique variants: ", uniqueN(coding[, .(CHROM, POS, REF, ALT)]),
        " | unique genes: ", uniqueN(coding$gene))

# --- cancer-gene overlap (any annotation, coding impact) ---
cg <- merge(long, cancer, by = "gene_key", allow.cartesian = TRUE)
cg_coding <- cg[is_coding == TRUE]
# prioritize: HCC driver > LoF HIGH > Act; order by is_hcc desc, impact, role
cg_coding[, impact_rank := fifelse(impact == "HIGH", 1L, 2L)]
cg_coding[, role_rank := fifelse(role == "LoF", 1L, fifelse(role == "Act", 2L, 3L))]
setorder(cg_coding, -is_hcc, impact_rank, role_rank, human_symbol)
fwrite(cg_coding[, .(human_symbol, rat_symbol, role, is_hcc, n_ctypes, ctypes,
                     CHROM, POS, REF, ALT, impact, annotation, cons_group,
                     hgvs_c, hgvs_p, feature_id, GT, DP, AD_ref, AD_alt, AF)],
       file.path(WD, "variants/cancer_coding_mutations.tsv"), sep = "\t")
message("[3] cancer-gene CODING mutations: ", nrow(cg_coding),
        " | unique cancer genes hit: ", uniqueN(cg_coding$human_symbol),
        " | HCC drivers hit: ", uniqueN(cg_coding[is_hcc == TRUE]$human_symbol))

# --- any-impact cancer overlap (for context: includes MODIFIER/LOW near cancer genes) ---
message("[4] cancer-gene ANY-impact annotations: ", nrow(cg),
        " | unique cancer genes: ", uniqueN(cg$human_symbol))

# --- Tp53 check ---
tp <- long[gene_key == "tp53" & is_coding == TRUE]
message("[5] Tp53 coding annotations recovered: ", nrow(tp))
if (nrow(tp)) print(tp[, .(CHROM, POS, REF, ALT, annotation, impact, hgvs_c, hgvs_p, AF)])

# --- consequence summary (primary = highest-impact annotation per variant) ---
prim <- long[order(factor(impact, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))),
             .SD[1], by = var_idx]
cons_tab <- prim[, .N, by = .(cons_group, impact)][order(-N)]
fwrite(cons_tab, file.path(WD, "variants/consequence_summary.tsv"), sep = "\t")
message("[6] consequence summary (primary annotation):")
print(cons_tab)

message("DONE")
