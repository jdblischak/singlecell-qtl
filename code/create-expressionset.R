#!/usr/bin/env Rscript

# Encapsulate information for each C1 chip into an ExpressionSet object.
#
# https://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf
#

# Setup ------------------------------------------------------------------------

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))

# Input arguments
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6, file.exists(args[1:5]))
f_molecules <- args[1]
f_lab <- args[2]
f_totals <- args[3]
f_verify <- args[4]
f_saf <- args[5]
f_rds <- args[6]

# Assay data: molecules counts -------------------------------------------------

molecules <- fread(sprintf("zcat %s", f_molecules), data.table = FALSE)

assay <- molecules %>% select(-sample:-well) %>% t
stopifnot(class(assay) == "matrix")
colnames(assay) <- molecules$sample
assay <- assay[order(rownames(assay)), ]

# Pheno data: single-cell metrics ----------------------------------------------

# Lab info
lab <- fread(f_lab, data.table = FALSE)
# Total counts
totals <- fread(f_totals, data.table = FALSE)
# Results from verifyBamID
verify <- fread(f_verify, data.table = FALSE)
stopifnot(lab$sample == totals$sample,
          lab$sample == verify$sample,
          lab$sample == colnames(assay))
# Calculate genes detected
detect_ce <- colSums(assay[grepl("WBGene", rownames(assay)), ] > 0)
detect_dm <- colSums(assay[grepl("FBgn", rownames(assay)), ] > 0)
detect_ercc <- colSums(assay[grepl("ERCC", rownames(assay)), ] > 0)
detect_hs <- colSums(assay[grepl("ENSG", rownames(assay)), ] > 0)

pheno <- cbind(lab %>% select(-sample),
               totals %>% select(-(sample:well)),
               detect_ce, detect_dm, detect_ercc, detect_hs,
               verify %>% select(-sample))
# Determine if predicted individual is one of the ones added to that C1 chip
pheno <- pheno %>% mutate(valid_id = chip_id %in% c(individual.1, individual.2,
                                                    individual.3, individual.4))
rownames(pheno) <- colnames(assay)

metadata <- data.frame(labelDescription = c(
  # Lab information
  "ID of C1 chip (i.e. processing date in MMDDYYYY)",
  "Well of C1 chip (96 total, rows A-H, cols 1-12)",
  "Batch the C1 chip was processed in (b1, b2, ...)",
  "The number of cells observed in the well via microscopy",
  "The cDNA concentration of the well prior to library prep",
  "Did the cell stain positive for TRA-1-60? (test of pluripotency)",
  "Individual # 1 included on this C1 chip",
  "Individual # 2 included on this C1 chip",
  "Individual # 3 included on this C1 chip",
  "Individual # 4 included on this C1 chip",
  "The concentration of spike-in added from D. Melanogaster",
  "The concentration of spike-in added from C. elegans",
  "The dilution factor of the ERCC spike-ins",
  "The set of indexes used for library prep (of the 3 sets of 96)",
  # Total counts
  "The number of raw reads",
  "The number of reads with a valid UMI",
  "The number of reads with a valid UMI that mapped to a genome",
  "The number of reads with a valid UMI that did *not* map to a genome",
  "The number of reads that mapped to the C. elegans genome",
  "The number of reads that mapped to the D. melanogaster genome",
  "The number of reads that mapped to the ERCC spike-in transcripts",
  "The number of reads that mapped to the H. sapiens genome",
  "The number of molecules (i.e. post UMI-deduplication)",
  "The number of molecules that mapped to the C. elegans genome",
  "The number of molecules that mapped to the D. melanogaster genome",
  "The number of molecules that mapped to the ERCC spike-in transcripts",
  "The number of molecules that mapped to the H. sapiens genome",
  # Number of genes detected
  "The number of C. elegans genes with at least one molecule",
  "The number of D. melanogaster genes with at least one molecule",
  "The number of ERCC genes with at least one molecule",
  "The number of H. sapiens genes with at least one molecule",
  # verifyBamID
  "verifyBamID: The predicted individual based on the sequencing data",
  "verifyBamID: chipmix is a metric for detecting sample swaps",
  "verifyBamID: freemix is a measure of contamination. 0 == good & 0.5 == bad",
  "verifyBamID: The number of SNPs that passed thresholds for AF and missingness",
  "verifyBamID: The number of sequences that overlapped SNPs",
  "verifyBamID: The average sequencing depth that covered a SNP",
  "verifyBamID: A minimun depth threshold for QC only (affects snps_w_min)",
  "verifyBamID: The number of SNPs that had the minimum depth (min_dp); QC only",
  "verifyBamID: Is the predicted individual 1 of the 4 added to the C1 chip?"
))
rownames(metadata) <- c(colnames(pheno))
pheno_anno <- new("AnnotatedDataFrame", data = pheno, varMetadata = metadata)
# Add back leading zero for experiment ID to get MMDDYYYY
pheno_anno$experiment <- sprintf("%08d", pheno_anno$experiment)
# concentration column must be numeric, even if all NA, so that it can be
# combined with other ExpressionSet objects
pheno_anno$concentration <- as.numeric(pheno_anno$concentration)

# Access info:
# sampleNames(pheno_anno)
# pData(pheno_anno)
# varMetadata(pheno_anno)

# Feature annotations ----------------------------------------------------------

exons <- fread(f_saf, data.table = FALSE)
colnames(exons) <- tolower(colnames(exons))
feature <- exons %>%
  group_by(geneid, chr, strand, name) %>%
  summarize(start = min(start), end = max(end)) %>%
  ungroup() %>%
  mutate(source = str_sub(chr, 1, 2),
         source = factor(source,
                         levels = c("ce", "dm", "ER", "hs"),
                         labels = c("C. elegans", "D. melanogaster",
                                    "ERCC", "H. sapiens")),
         source = as.factor(source)) %>%
  arrange(geneid)
geneid <- feature$geneid
feature <- feature %>%
  select(chr, start, end, name, strand, source) %>%
  as.data.frame()
rownames(feature) <- geneid
feature_meta <- data.frame(labelDescription = c(
  "Chromosome",
  "Most 5' start position (GRCh37/hg19; 1-based; inclusive)",
  "Most 3' end position (GRCh37/hg19; 1-based; inclusive)",
  "Gene name",
  "Strand (+ = positive/forward; - = negative/reverse)",
  "Source of RNA"
))
rownames(feature_meta) <- colnames(feature)

feature_anno <- new("AnnotatedDataFrame", data = feature,
                    varMetadata = feature_meta)

# Experiment data --------------------------------------------------------------

experiment <- new(
  "MIAME",
  name = "Abhishek K. Sarkar, Po-Yuan Tung, John D. Blischak, Yang I. Li, Matthew Stephens, and Yoav Gilad",
  lab = "Gilad (University of Chicago)",
  contact = "https://github.com/jdblischak/singlecell-qtl/issues",
  title = "Discovery and characterization of variance QTLs in human induced pluripotent stem cells",
  abstract =
"Quantification of gene expression levels at the single cell level has
revealed that gene expression can vary substantially even across a
population of homogeneous cells. However, it is currently unclear what
genomic features control variation in gene expression levels, and whether
common genetic variants may impact gene expression variation. Here, we
take a genome-wide approach to identify expression variance quantitative
trait loci (vQTLs). To this end, we generated single cell RNA-seq
(scRNA-seq) data from induced pluripotent stem cells (iPSCs) derived from
53 Yoruba individuals. We collected data for a median of 95 cells per
individual and a total of 5,447 single cells, and identified 241 mean
expression QTLs (eQTLs) at 10% FDR, of which 82% replicate in bulk RNA-seq
data from the same individuals. We further identified 14 vQTLs at 10% FDR,
but demonstrate that these can also be explained as effects on mean
expression. Our study suggests that dispersion QTLs (dQTLs), which could
alter the variance of expression independently of the mean, have
systematically smaller effect sizes than eQTLs. We estimate that at least
300 cells per individual and 400 individuals would be required to have
modest power to detect the strongest dQTLs in iPSCs. These results will
guide the design of future studies on understanding the genetic control of
gene expression variance.",
  url = "https://jdblischak.github.io/singlecell-qtl",
  preprocessing = list(
    umi_extract =
"We used `umi_tools extract` (UMI-tools 0.5.3) to extract the 6 base
pair unique molecular identifier (UMI) from the 5’ end of each read.",
    mapping =
"We used `subjunc` (Subread 1.5.3) to map reads to the human genome
(Ensembl GRCh37.75; chromosomes 1-22, X, Y, MT) and ERCC spike-ins
(http://tools.invitrogen.com/downloads/ERCC92.fa). Some samples included
spike-in RNA from Drosophila melanogaster or Caenorhabditis elegans, so
we also included the Ensembl genomes BDGP5.75 and WBcel235.75.",
    counting =
"We used `featureCounts` (Subread 1.5.3) to count the number of reads
for all protein-coding genes (Ensembl GRCh37 release 75) and the ERCC
spike-in genes. We performed strand-specific counting (flag -s 1)
because the UMI protocol preserves sequence strand information.",
    umi_dedup =
"We used `umi_tools dedup` (UMI-tools 0.5.3) to deduplicate reads with
the same UMI and start position to molecules.",
    verify =
"We used `verifyBamID` (1.1.3) to identify the individual of origin for
each sample based on the overlap of the RNA-seq reads with the known
genotypes."
  ))

# Create ExpressionSet object --------------------------------------------------

eset <- ExpressionSet(assayData = assay,
                      phenoData = pheno_anno,
                      featureData = feature_anno,
                      experimentData = experiment)

saveRDS(eset, file = f_rds)

# Access info:
# eset
# dim(eset)
# sampleNames(eset)
# phenoData(eset)
# varMetadata(eset)
# featureData(eset)
# fvarMetadata(eset)
# experimentData(eset)
# abstract(eset)
# preproc(eset)
