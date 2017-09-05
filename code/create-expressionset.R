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

pheno <- cbind(lab %>% select(-sample),
               totals %>% select(-(sample:well)),
               verify %>% select(-sample))
# Determine if predicted individual is one of the ones added to that C1 chip
pheno <- pheno %>% mutate(valid_id = chip_id %in% c(individual.1, individual.2,
                                                    individual.3, individual.4))
rownames(pheno) <- colnames(assay)

metadata <- data.frame(labelDescription = c(
  # Lab information
  "ID of C1 chip (i.e. processing date in MMDDYYYY)",
  "Well of C1 chip (96 total, rows A-H, cols 1-12)",
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
# Add back leading zero for experiment ID
pheno_anno$experiment <- paste0("0", pheno_anno$experiment)
# ERCC column must be character, even if all NA, so that it can be combined with
# other ExpressionSet objects
pheno_anno$ERCC <- as.character(pheno_anno$ERCC)
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
  name = "PoYuan Tung, Harold Pimentel, Joyce Hsiao, John Blischak",
  lab = "Gilad/Pritchard",
  contact = "https://github.com/jdblischak/singlecell-qtl/issues",
  title = "singlecell-qtl",
  abstract = "Mapping genetic variants associated with cell-to-cell variation in gene expression in the HapMap Yoruba population.",
  url = "https://jdblischak.github.io/singlecell-qtl",
  preprocessing = list(
    mapping = "Mapped to GRCh37 with Sunjunc",
    umi = "UMIs were extracted and deduplicated with umi_tools",
    counting = "Counts per gene (Ensembl 75) were obtained with featureCounts",
    verify = "Individuals were identified by verifyBamID"
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
