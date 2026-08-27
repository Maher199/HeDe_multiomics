#!/usr/bin/env Rscript
# fig7: genome-browser views of top-5 driver genes (coverage + gene model + domains + lollipops)
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(svglite))

WD <- "/workspace/hede_followup"
FD <- file.path(WD, "fig7_data")
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

genes <- c("TP53","KMT2D","EP300","FAT1","DICER1")
ginfo <- list(
  TP53  =list(chr="chr10", start=54798871,  end=54810300,  strand="+"),
  KMT2D =list(chr="chr7",  start=131859696, end=131901032, strand="-"),
  EP300 =list(chr="chr7",  start=114987857, end=115058652, strand="+"),
  FAT1  =list(chr="chr17", start=53909759,  end=54029175,  strand="-"),
  DICER1=list(chr="chr6",  start=129392298, end=129457252, strand="-")
)
# key domain labels to show per gene (regex on description)
key_dom <- list(
  TP53  ="DNA binding",
  KMT2D ="^SET$",
  EP300 ="HAT|Bromo|KIX",
  FAT1  ="Cadherin 20$",
  DICER1="RNase III|PAZ"
)
# mutations to label per gene (hgvs_p); KMT2D limited to avoid clutter
label_mut <- list(
  TP53  ="p.Arg271Ser",
  KMT2D =c("p.Thr4847fs","p.Gly1866Ser","p.Gly1933Ser"),
  EP300 ="p.Gln313*",
  FAT1  ="p.Thr2257Ile",
  DICER1="p.Glu1701Lys"
)

mut <- fread(file.path(FD, "mutations.tsv"))

build_gene <- function(g) {
  gi <- ginfo[[g]]
  x0 <- gi$start/1e6; x1 <- gi$end/1e6
  pad <- (x1-x0)*0.02
  xl <- c(x0-pad, x1+pad)

  cov <- fread(file.path(FD, paste0("cov_", g, ".tsv"))); cov[, mb:=pos/1e6]
  ex  <- fread(file.path(FD, paste0("exons_", g, ".tsv")))
  cd  <- fread(file.path(FD, paste0("cds_", g, ".tsv")))
  dom <- fread(file.path(FD, paste0("domains_", g, ".tsv")))
  dom <- dom[ftype %in% c("Domain","Zinc finger","DNA binding")]
  dom[, `:=`(g0=gstart/1e6, g1=gend/1e6)]
  dom[, dlabel := fifelse(is.na(description) | description=="", ftype, description)]
  dom[, key := grepl(key_dom[[g]], dlabel)]
  # short display labels for key domains
  dom[, short := dlabel]
  dom[grepl("DNA binding", dlabel), short := "DNA-binding"]
  dom[grepl("CBP/p300-type HAT", dlabel), short := "HAT"]
  dom[grepl("Cadherin 20", dlabel), short := "Cadherin 20"]
  mm <- mut[human_symbol==g]; mm[, mb:=POS/1e6]
  mm[, lab := ifelse(hgvs_p %in% label_mut[[g]], gsub("p\\.","",hgvs_p), "")]
  # dodge labels: stagger close lollipops vertically to avoid overlap
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
    labs(y="AF", title=paste0(g, "  (", gi$chr, ", ", gi$strand, " strand)")) +
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

  # 4) domains (label only key ones, placed above boxes to avoid cramping)
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

w <- 8.5; h <- 13
svglite(file.path(FIGD, "fig7_driver_gene_browser.svg"), width=w, height=h)
print(fig); dev.off()
ggsave(file.path(FIGD, "fig7_driver_gene_browser.png"), fig, width=w, height=h, dpi=150)
cat("saved fig7 (", w, "x", h, ")\n")
