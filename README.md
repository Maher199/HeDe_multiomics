# HeDe PacBio HiFi multi-omic dataset (GRCr8)

Processed data accompanying the Data Descriptor: *"A PacBio HiFi multi-omic dataset of
the HeDe rat hepatocellular carcinoma cell line"* (Scientific Data, under review).

This archive contains the derived data products from deep PacBio HiFi sequencing of the
HeDe rat hepatocellular carcinoma (HCC) cell line: a de novo assembly, small-variant and
structural-variant calls, copy-number segments, and genome-wide 5-methylcytosine (5mC)
profiles, plus the analysis code. Raw HiFi reads are deposited separately in the NCBI
Sequence Read Archive (SRA) under BioProject accession [![NCBI SRA](https://img.shields.io/badge/NCBI%20SRA-PRJNA1517565-blue.svg)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1517565/).

**Reference genome:** all coordinates use the rat reference **GRCr8**
(GCF_036323735.1; UCSC rn8). Chromosomes are RefSeq `NC_*` accessions
(NC_086019.1 = chr1 … NC_086038.1 = chr20; NC_086039.1 = chrX; NC_086040.1 = chrY).

**Sequencing:** PacBio Revio, circular consensus sequencing (CCS 8.2.0); base
modifications (5mC and 6mA) called with jasmine 2.3.0; demultiplexing with lima 2.11.0.

**Version 2 addition (2026-08):** the `wgs/` directory adds the derived products of the
companion Illumina whole-genome sequencing (WGS) dataset — one HeDe tumor library (H_1)
and one strain-matched healthy F344 liver library (RL_1, different donor animal) —
including somatic and germline variant calls, copy-number segments, cross-platform
concordance tables, and mutational-signature refits. Raw WGS reads are available under BioProject accession
`PRJNA1517565`.

---

## Directory contents

### assembly/
- `HeDe_assembly.fasta.gz` — de novo assembly (hifiasm), 605 primary contigs,
  2.81 Gb total, contig N50 65.2 Mb. Contig names use hifiasm `ptgNNNNNNl` convention.

### variants/
- `HeDe_novel_PASS.vcf.gz` (+ `.tbi`) — 3,319,205 novel PASS small variants (SNVs and
  indels) called with Clair3 2.0.2 and filtered to PASS; "novel" = absent from rat dbSNP.
- `HeDe_novel_PASS.ann.vcf.gz` — the same variants annotated with SnpEff 5.4 (custom
  GRCr8 database); functional consequences in the `ANN` INFO field.

### sv/
- `HeDe_sniffles2.vcf.gz` (+ `.tbi`) — structural variants called with Sniffles2 2.8.0
  (`--minsvlen 50 --output-rnames`). 73,273 VCF records in total, of which 72,831 are PASS
  (DEL 38,399; INS 32,599; BND 1,308; INV 295; DUP 230). Analyses in the accompanying
  manuscripts use the 72,831 PASS calls.
- `HeDe_sv_table.tsv` — flattened table of the 72,831 PASS calls. Tab-separated, **no header
  row**, four columns: `chrom`, `pos`, `svtype`, `svlen` (svlen is negative for deletions).

### cnv/
- `HeDe_cnv_20kb_bins.tsv` — read-depth in 20-kb bins (mosdepth 0.3.14).
- `HeDe_cnv_segments.tsv` — copy-number segments from circular binary segmentation
  (DNAcopy 1.80.0). Columns: `ID, chr, start, end, nmark, log2r, width, state`
  (`log2r` = log2 copy-number ratio; `state` = gain/loss/neutral).

### methylation/
- `HeDe_cpg_5mC.bed.gz` — per-CpG 5mC from modkit 0.6.4 (bedMethyl). Key columns:
  chrom, start, end, mod code (`m`), read coverage (col 5), percent methylated (col 11).
- `HeDe_promoter_5mC.tsv` — coverage-weighted mean 5mC per promoter (2 kb upstream).
- `cgi_5mC.tsv` — mean 5mC per CpG island. 15,855 rows: of the 16,235 CGIs annotated in GRCr8
  by Gardiner-Garden criteria, these are the ones with callable coverage (>= 10x). The
  `meth_state` column uses <= 0.2 hypomethylated / >= 0.6 hypermethylated, giving
  11,513 hypomethylated / 803 intermediate / 3,539 hypermethylated.
- `retrotransposon_5mC.tsv` — mean 5mC per retrotransposon family (LINE/SINE/LTR).
- `imprinted_5mC.tsv` — mean 5mC at 29 imprinted control regions.

### tables/
- `coding_mutations_all.tsv` — all 9,508 unique coding variants (HIGH + MODERATE impact).
- `cancer_coding_mutations.tsv` — the 72 unique coding variants in 55 cancer genes.
- `cancer_mutations_functional_impact.tsv` — per-variant functional annotation: gene,
  driver role, intOGen HCC-driver status, GRCr8 coordinate, consequence, protein change,
  genotype/allele fraction, rat SIFT, human-ortholog change with CADD and AlphaMissense
  scores, overlapping CNV state, LOH class, and a composite functional call.
- `promoter_5mC_x_cnv_cancer.tsv` — promoter 5mC integrated with copy-number state for
  cancer genes.
- `mutational_spectrum_6class.csv`, `mutational_spectrum_96.csv` — SNV spectra.
- `signature_contributions_targeted.csv`, `HeDe_signature_contributions.csv` — COSMIC
  mutational-signature contributions.
- `HeDe_karyotype_by_chromosome.tsv` — per-chromosome molecular karyotype summarising the
  DNAcopy segments in `cnv/HeDe_cnv_segments.tsv`: segment count, gain/loss/neutral base
  pairs, width-weighted mean log2 ratio, chromosome length, covered length, the three
  fractions, and a whole-chromosome call. This is the table reported in the accompanying
  manuscripts (chr13 gain, chr12 and chrY loss, chr2 partial gain, chr19/chr16/chr5 partial
  loss, the other 15 chromosomes neutral). Fractions are of total chromosome length, not of
  the covered length. Calls use 0.50 for whole-chromosome and 0.15 for partial; note that
  any partial cut-off above 0.074 and up to 0.172 gives the identical set of calls, because
  no chromosome has an altered fraction in that interval (highest neutral chromosome is
  chr18 at 0.074, lowest partial chromosome is chr5 at 0.172).
- `GRCr8_chromosome_lengths.tsv` — the 22 assembled GRCr8 chromosome lengths with their
  RefSeq accessions (`chr`, `accession`, `length`), taken from the reference FASTA index.
  Used by `code/make_karyotype.R` so it runs from this archive without the reference.
- `read_lengths_full.txt.gz` — one integer read length per line, in BAM order, for all
  5,565,830 HiFi reads (gzip, 13 MB; 33 MB uncompressed). Recovered from the raw BAM index
  `H_mix.hifi_reads.bam.pbi`, so the read-length statistics can be checked without
  downloading the 76 GB BAM.
- `read_length_summary.tsv` — yield and length statistics over the full read set: total
  128,127,172,081 bp (128.13 Gb), mean 23,020.3 bp, median 22,601 bp, N50 23,294 bp,
  N90 18,856 bp, range 108-62,043 bp, plus the 1/5/25/75/95/99th percentiles. These are the
  values reported in the accompanying manuscripts. Note that a head-of-file subsample is
  length-biased and will not reproduce them (the first 500,000 reads give an N50 of
  23,526 bp), so always use the complete distribution.

### code/
Analysis scripts (R and Python) used to generate the tables and figures.
`read_lengths_from_pbi.R` regenerates the two read-length files above directly from the
raw BAM index: `Rscript code/read_lengths_from_pbi.R <in.bam.pbi> tables/`.
`make_karyotype.R` regenerates `tables/HeDe_karyotype_by_chromosome.tsv` and the molecular
karyotype ideogram from `cnv/HeDe_cnv_segments.tsv` and `tables/GRCr8_chromosome_lengths.tsv`,
both of which are in this archive, so it runs standalone: `Rscript code/make_karyotype.R`.
This script was reconstructed after its working copy was lost; it is verified by output
rather than byte-recovered, and reproduces the deposited karyotype table exactly (MD5
`10a052164fa27895eefa5d6a27755050`). Its header documents the one cosmetic parameter and
the one threshold that the surviving outputs do not fully constrain.

### wgs/
Derived products of the companion Illumina WGS dataset: HeDe tumor DNA (sample H_1,
NucleoSpin extraction, 90.1 Gb) and a strain-matched healthy F344 liver from a different
donor animal (sample RL_1, phenol-chloroform extraction, 101.6 Gb). Novogene NovaSeq,
paired-end 150 bp, PCR-amplified libraries. Reads were aligned to GRCr8 with minimap2
2.28 (`-ax sr`) and duplicate-marked with samtools. H_1: 98.15% of reads mapped, 32.6x
raw / 19.3x deduplicated effective coverage. RL_1: 96.07% mapped, 38.6x raw / 20.7x
effective. Because the liver control comes from a different F344 donor, germline variants
private to the HeDe donor animal cannot be distinguished from somatic variants; the raw
somatic call set is therefore an upper bound, as disclosed in the accompanying manuscripts.

- `HeDe_WGS_somatic_mutect2.vcf.gz` (+ `.tbi`) — somatic SNV/indel calls from Mutect2
  (GATK 4.6.1.0), tumor H_1 vs liver RL_1, after FilterMutectCalls: 537,751 records,
  82,545 PASS (51,374 SNVs and 31,171 indels). No germline resource or panel of normals
  was used (none with allele frequencies exists for rat), so the PASS set includes
  donor-private germline variants (see above). The HeDe driver mutation Tp53
  p.Arg271Ser (NC_086028.1:54,807,961 C>A) is recovered at PASS, tumor AF 0.397.
- `HeDe_WGS_liver_germline.vcf.gz` (+ `.tbi`) — germline variants of the RL_1 liver
  animal against GRCr8 (FreeBayes 1.0.2): 8,236,556 records, of which 4,504,427 are PASS
  (the remainder are LOWQUAL and retained for transparency). An F344-vs-reference
  germline catalog.
- `HeDe_WGS_manta_somaticSV.vcf.gz` (+ `.tbi`) — 2,329 somatic SV records (Manta 1.6.0,
  tumor vs liver), of which 44 are PASS.
- `HeDe_WGS_H1_20kb_bins.bed.gz`, `HeDe_WGS_RL1_20kb_bins.bed.gz` — per-sample read depth
  in 20-kb bins (mosdepth); input to the WGS copy-number analysis.
- `HeDe_WGS_cnv_segments.tsv` — DNAcopy circular-binary-segmentation output on the
  tumor/liver log2 ratio (same method and binning as the PacBio `cnv/` analysis).
- `HeDe_WGS_karyotype_by_chromosome.tsv` — per-chromosome summary of those segments with
  the same 0.50 whole-chromosome / 0.15 partial call thresholds used for the PacBio
  karyotype. At these sensitive thresholds many chromosomes pick up minor partial calls;
  the authoritative cross-platform comparison is the next file.
- `HeDe_karyotype_pacbio_vs_wgs.tsv` — per-chromosome comparison of the PacBio
  width-weighted segment mean vs the robust WGS median-20kb-bin log2(tumor/liver)
  estimator. Concordant: chr13 gain, chrY loss, chr5 and chr16 partial losses, and 14
  neutral chromosomes. Discordant and resolved by the adjudication below: chr12, chr10,
  chr19 (and the weak PacBio chr2 partial-gain call, wlog2r +0.133, is not supported by
  WGS, -0.068).
- `HeDe_chr10_12_19_verdicts.tsv` — final per-chromosome copy-number verdicts with all
  supporting metrics. chr12, chr10 and chr19 are diploid-neutral in the tumor. The
  PacBio chr12-loss call is a long-read mappability artifact (chr12 is the rDNA/NOR
  chromosome; its duplicative content is tandemly arrayed, not segmental), and the WGS
  single-sample window-mean inflations are the opposite-sign short-read pile-up
  artifacts; the tumor/normal median-bin ratio is the robust estimator.
- `adjudication_windows_chr10_12_19.tsv`, `adjudication_chr_summary.tsv`,
  `window_profiles_chr10_12_19.tsv`, `mapq_profiles_500kb.tsv` — the per-window evidence
  behind the verdicts (WGS tumor and liver depth, tumor/normal median-bin log2 ratio,
  PacBio relative depth, tandem-repeat fractions, MAPQ profiles).
- `segdup_selfaln_summary.tsv` — de novo segmental-duplication screen (minimap2
  self-alignment of chr9/chr10/chr12/chr19): diagonal-only alignments for all four
  chromosomes, i.e. no large interspersed segmental duplications; duplicative content is
  tandemly arrayed. GRCr8 carries no curated segmental-duplication annotation.
- `HeDe_WGS_somatic.96spectrum.csv`, `HeDe_WGS_somatic.6class.csv` — 96-channel
  trinucleotide and 6-class spectra of the 51,374 PASS SNVs (100% REF-vs-FASTA
  validated).
- `HeDe_WGS_somatic.AFlt35.96spectrum.csv`, `HeDe_WGS_somatic.AFge35.96spectrum.csv` —
  the same spectra split at tumor allele fraction 0.35 (subclonal-enriched vs
  clonal-enriched).
- `HeDe_WGS_signature_targeted.csv`, `HeDe_WGS_signature_contributions.csv` — COSMIC
  v3.3 NNLS refits (14-signature targeted panel and unrestricted).
- `signature_fit_pacbio_vs_wgs.csv`, `signature_fit_AFstrata_vs_pacbio.tsv` —
  cross-platform signature comparison. The alkylation signature SBS11 (DMN is an
  alkylating nitrosamine) contributes 37.9% of the clonal (AF >= 0.35) WGS spectrum but
  only 3.6% of the tumor-only PacBio novel-variant fit, where it is diluted by the
  germline/background variant mass.
- `sv_concordance_by_type.tsv`, `sv_concordance_by_size.tsv`,
  `sv_concordance_merged.tsv` — Manta vs Sniffles2 concordance (SURVIVOR, type-matched,
  +/-500 bp). Only 2 of the 44 Manta PASS SVs have a Sniffles2 partner (both large
  deletions); the two SV sets are otherwise technology-specific.
- `HeDe_WGS_at_pacbio_sites.vcf.gz` (+ `.tbi`) — WGS genotypes at 973,944 PacBio
  novel-SNV sites (targeted genotyping of the tumor WGS).
- `HeDe_WGS_vs_pacbio_concordance_by_depth.tsv`,
  `HeDe_WGS_vs_pacbio_concordance_merged.tsv.gz` — cross-platform genotype concordance
  stratified by WGS depth (PacBio ALT detected in 69.9% of sites at >= 20x WGS depth).
- `fastqc_gate_summary.tsv` — FastQC module status for all WGS lanes of both samples.
- `sample_manifest.tsv` — sample and tube metadata for the two WGS libraries.

---

## Checksums
`MD5SUMS` lists MD5 checksums for all files; `MANIFEST.tsv` lists file sizes.

## Citation
If you use these data, please cite the Data Descriptor (see above) and this dataset
[![DOI](https://zenodo.org/badge/854581223.svg)](https://doi.org/10.5281/zenodo.22126105)
