#!/usr/bin/env Rscript
# fig8: genome-browser views of second-tier intOGen HCC driver genes
# (coverage + gene model + domains + lollipops)
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(svglite))

WD <- "/workspace/hede_followup"
FD <- file.path(WD, "fig8_data")
FIGD <- file.path(WD, "figures")
dir.create(FIGD, showWarnings = FALSE)

# Okabe-Ito
OI <- c(orange="#E69F00", sky="#56B4E9", green="#009E73", yellow="#F0E442",
        blue="#0072B2", vermillion="#D55E00", purple="#CC79A7", black="#000000")
impact_col <- c("HIGH"=OI[["vermillion"]], "MODERATE"=OI[["orange"]])

base_theme <- theme_bw(base_family="Liberation Sans") +
  theme(text=element_text(family="Liberation Sans"),
        panel.grid=element_blank(),
        plot.margin=margin(2,8,2,8))

# ordered by intOGen n_ctypes (driver strength)
genes <- c("KMT2C","ARID1B","TET1","CLTC","WNK2","DHX9")
ginfo <- list(
  KMT2C =list(chr="chr4",  start=10353698, end=10755965, strand="+", nct=32),
  ARID1B=list(chr="chr1",  start=47973199, end=48328793, strand="+", nct=11),
  TET1  =list(chr="chr20", start=25766806, end=25839598, strand="+", nct=5),
  CLTC  =list(chr="chr10", start=72014984, end=72073308, strand="-", nct=3),
  WNK2  =list(chr="chr17", start=15709648, end=15818874, strand="+", nct=3),
  DHX9  =list(chr="chr13", start=68152813, end=68189580, strand="-", nct=1)
)
# extra domain-track rows to include beyond Domain/Zinc finger/DNA binding (regex on description)
include_dom <- list(
  KMT2C="", ARID1B="", TET1="", CLTC="CHCR|Trimerization", WNK2="", DHX9="^RGG"
)
# key domain(s) to highlight + label (regex on description)
key_dom <- list(
  KMT2C="^SET$", ARID1B="^ARID$", TET1="CXXC",
  CLTC="Trimerization", WNK2="Protein kinase", DHX9="^Helicase|^RGG"
)
# mutations to label (all; 1-2 per gene)
label_mut <- list(
  KMT2C=c("p.Gln1808del","p.Gln3310del"),
  ARID1B=c("p.Ala1037Thr","p.Gly1948Asp"),
  TET1="p.Val926Ile",
  CLTC="p.Ser1650Asn",
  WNK2="p.Glu177dup",
  DHX9="p.Gly1390dup"
)

mut <- fread(file.path(FD, "mutations.tsv"))

build_gene <- function(g) {
  gi <- ginfo[[g]]
  x0 <- gi$start/1e6; x1 <- gi$end/1e6
  pad <- (x1-x0)*0.05
  xl <- c(x0-pad, x1+pad)

  cov <- fread(file.path(FD, paste0("cov_", g, ".tsv"))); cov[, mb:=pos/1e6]
  ex  <- fread(file.path(FD, paste0("exons_", g, ".tsv")))
  cd  <- fread(file.path(FD, paste0("cds_", g, ".tsv")))
  dom <- fread(file.path(FD, paste0("domains_", g, ".tsv")))
  inc <- include_dom[[g]]
  keep <- dom$ftype %in% c("Domain","Zinc finger","DNA binding")
  if (nzchar(inc)) keep <- keep | grepl(inc, dom$description)
  dom <- dom[keep]
  dom[, `:=`(g0=gstart/1e6, g1=gend/1e6)]
  dom[, dlabel := fifelse(is.na(description) | description=="", ftype, description)]
  dom[, key := grepl(key_dom[[g]], dlabel)]
  dom[, short := dlabel]
  dom[grepl("Protein kinase", dlabel), short := "Kinase"]
  dom[grepl("Helicase ATP-binding", dlabel), short := "Helicase-ATP"]
  dom[grepl("Helicase C-terminal", dlabel), short := "Helicase-C"]
  mm <- mut[human_symbol==g]; mm[, mb:=POS/1e6]
  mm[, lab := ifelse(hgvs_p %in% label_mut[[g]], gsub("p\\.","",hgvs_p), "")]
  mm[lab!="", lab_y := AF + 0.10]
  labd <- mm[lab!=""][order(mb)]
  if (nrow(labd) > 1) {
    for (i in 2:nrow(labd)) {
      if ((labd$mb[i]-labd$mb[i-1]) < (x1-x0)*0.10) labd$lab_y[i] <- labd$lab_y[i-1] + 0.20
    }
    mm[lab!="", lab_y := labd$lab_y[match(hgvs_p[lab!=""], labd$hgvs_p)]]
  }

  # 1) lollipop
  p_lol <- ggplot(mm) +
    geom_segment(aes(x=mb, xend=mb, y=0, yend=AF, color=impact), linewidth=0.5) +
    geom_point(aes(x=mb, y=AF, color=impact), size=2.6) +
    geom_text(data=mm[lab!=""], aes(x=mb, y=lab_y, label=lab),
              size=2.4, family="Liberation Sans", vjust=0) +
    scale_color_manual(values=impact_col, guide="none") +
    scale_y_continuous(limits=c(0,1.35), breaks=c(0,0.5,1), expand=expansion(0)) +
    coord_cartesian(xlim=xl) +
    labs(y="AF", title=paste0(g, "  (", gi$chr, ", ", gi$strand, " strand; intOGen ", gi$nct, " ctypes)")) +
    base_theme +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=7), axis.title.y=element_text(size=8),
          plot.title=element_text(size=10, face="bold", margin=margin(b=1)))

  # 2) coverage
  p_cov <- ggplot(cov, aes(x=mb, y=depth)) +
    geom_area(fill=OI[["sky"]], alpha=0.6, linewidth=0.2, color=OI[["blue"]]) +
    scale_y_continuous(breaks=scales::pretty_breaks(3), expand=expansion(c(0,0.05))) +
    coord_cartesian(xlim=xl) +
    labs(y="Cov") + base_theme +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=6.5), axis.title.y=element_text(size=8))

  # 3) gene model
  p_gene <- ggplot() +
    geom_segment(aes(x=x0, xend=x1, y=0, yend=0), linewidth=0.5, color="grey40") +
    geom_rect(data=ex, aes(xmin=start/1e6, xmax=end/1e6, ymin=-0.35, ymax=0.35),
              fill="grey70", color=NA) +
    geom_rect(data=cd, aes(xmin=start/1e6, xmax=end/1e6, ymin=-0.6, ymax=0.6),
              fill="grey20", color=NA) +
    annotate("text", x=x1, y=0.9, label=ifelse(gi$strand=="+","5'\u21923'","3'\u21905'"),
             size=2.6, family="Liberation Sans", hjust=1) +
    scale_y_continuous(limits=c(-1,1), expand=expansion(0)) +
    coord_cartesian(xlim=xl) +
    labs(y="Gene") + base_theme +
    theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(),
          panel.border=element_blank(), axis.line=element_blank())

  # 4) domains
  p_dom <- ggplot(dom) +
    geom_rect(aes(xmin=g0, xmax=g1, ymin=0, ymax=1, fill=key), color="white", linewidth=0.2) +
    geom_text(data=dom[key==TRUE], aes(x=(g0+g1)/2, y=1.12, label=short),
              size=2.1, family="Liberation Sans", color="black", vjust=0) +
    scale_fill_manual(values=c("TRUE"=OI[["green"]], "FALSE"="grey82"), guide="none") +
    scale_y_continuous(limits=c(0,1.5), expand=expansion(0)) +
    coord_cartesian(xlim=xl) +
    labs(y="Dom", x=paste0(gi$chr, " position (Mb)")) + base_theme +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          axis.title=element_text(size=8), axis.text.x=element_text(size=7))

  p_lol / p_cov / p_gene / p_dom + plot_layout(heights=c(2.1,1.1,0.7,0.7))
}

panels <- lapply(genes, build_gene)
fig <- wrap_plots(panels, ncol=1) +
  plot_annotation(theme=theme(text=element_text(family="Liberation Sans")))

w <- 8.5; h <- 15.5
svglite(file.path(FIGD, "fig8_secondtier_driver_browser.svg"), width=w, height=h)
print(fig); dev.off()
ggsave(file.path(FIGD, "fig8_secondtier_driver_browser.png"), fig, width=w, height=h, dpi=150)
cat("saved fig8 (", w, "x", h, ")\n")
