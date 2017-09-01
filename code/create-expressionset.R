#!/usr/bin/env Rscript

# Encapsulate information for each C1 chip into an ExpressionSet object.
#
# https://www.bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf
#

library("Biobase")
library("data.table")
library("dplyr")
library("stringr")

molecules <- fread("zcat ../data/molecules/03302017.txt.gz", data.table = FALSE)

# Assay data: molecules counts -------------------------------------------------

assay <- molecules %>% select(-sample:-well) %>% t
stopifnot(class(assay) == "matrix")
colnames(assay) <- molecules$sample
assay <- assay[order(rownames(assay)), ]

# Pheno data: single-cell metrics ----------------------------------------------

# Total counts
totals <- fread("../data/totals/03302017.txt", data.table = FALSE)
# Results from verifyBamID
verify <- fread("../data/verify/03302017.txt", data.table = FALSE)
stopifnot(totals$sample == verify$sample,
          totals$sample == colnames(assay))

pheno <- merge(totals, verify) %>% select(-sample)
rownames(pheno) <- colnames(assay)

metadata <- data.frame(labelDescription = c(colnames(pheno)))
rownames(metadata) <- c(colnames(pheno))
pheno_anno <- new("AnnotatedDataFrame", data = pheno, varMetadata = metadata)

# Access info:
# sampleNames(pheno_anno)
# pData(pheno_anno)
# varMetadata(pheno_anno)

# Feature annotations ----------------------------------------------------------

saf <- "/scratch/midway2/jdblischak/genome-ensembl-release-75/exons-ensembl-release-75.saf"
exons <- fread(saf, data.table = FALSE)
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

saveRDS(eset, "../data/eset/03302017.rds")
