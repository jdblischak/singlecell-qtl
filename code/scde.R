#!/usr/bin/env Rscript

# Submit comptutationally intensive SCDE/Pagoda steps to cluster.
#
# Creates the following three files:
#
# * ../data/eset-sub-knn{append}.rds - returned by knn.error.models
# * ../data/pwpca-eset-dep{append}.rds - returned by pagoda.pathway.wPCA
# * ../data/clpca-eset-dep{append}.rds - returned by pagoda.gene.clusters
#
# where {append} is set inside the script to designate the input data.
#
# Delete them first if you want to recreate them with this script.
#
# Submit from code/ directory:
#
# sbatch --job-name=scde --output=scde-slurm.out --time=24:00:00 --partition=bigmem2 --mem=128G --tasks-per-node=8 scde.R

library("biomaRt")
library("cowplot")
library("edgeR")
library("ggplot2")
library("knitr")
library("scde")
theme_set(theme_cowplot())
source("../code/functions.R")
library("Biobase") # has to be loaded last to use `combine`
library("GO.db")

append <- "-filtered-b2-b5"

# Setup ------------------------------------------------------------------------
cat(sprintf("\n\nSetup...\n\n"))

# Import data.
eset <- readRDS("../data/eset.rds")

# Keep human genes and ERCC
eset <- eset[fData(eset)$source %in% c("H. sapiens", "ERCC") , ]

# Only keep high-quality single cells.
quality <- read.table("../data/quality-single-cells.txt", stringsAsFactors = FALSE)
colnames(quality) <- c("sample", "quality")
eset <- eset[, quality$quality]

# Only keep genes that passed the filters
genes <- read.table("../data/genes-pass-filter.txt", stringsAsFactors = FALSE)
colnames(genes) <- c("gene", "passed")
eset <- eset[genes$passed, ]

eset_data <- exprs(eset)

# Limit cells to batches 2-5 (not all cells in batch 1 had ERCC spike-in added)
eset_data_sub <- eset_data[, pData(eset)$batch != "b1"]

# SCDE -------------------------------------------------------------------------
cat(sprintf("\n\nFit error models...\n\n"))

fname_error_model <- paste0("../data/eset-sub-knn", append, ".rds")

if (file.exists(fname_error_model)) {
  eset_sub_knn <- readRDS(fname_error_model)
} else {
  ## Fitting error models
  eset_sub_knn <- knn.error.models(eset_data_sub,
                                   k = round(ncol(eset_data_sub) / 4),
                                   n.cores = 8,
                                   min.count.threshold = 1,
                                   min.nonfailed = 5,
                                   save.model.plots = FALSE)

  saveRDS(eset_sub_knn, fname_error_model)
}

## Normalizing variance
varinfo_eset_sub <- pagoda.varnorm(eset_sub_knn,
                                   counts = eset_data_sub,
                                   trim = 3,
                                   max.adj.var = 5,
                                   n.cores = 8,
                                   plot = FALSE)

## Controlling for sequencing depth
varinfo_dep <- pagoda.subtract.aspect(varinfo_eset_sub,
                                      colSums(eset_data_sub[, rownames(eset_sub_knn)]>0))

# Prepare GO terms -------------------------------------------------------------
cat(sprintf("\n\nPrepare GO terms...\n\n"))

# Initialize the connection to the Ensembl BioMart Service
# Available datasets can be listed with
# listDatasets(useMart("ENSEMBL_MART_ENSEMBL", host = "feb2014.archive.ensembl.org"))
# Use mmusculus_gene_ensembl for mouse
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl",
                   host = "feb2014.archive.ensembl.org")

# Constructs a dataframe with two columns: hgnc_symbol and go_id
# If rownames are Ensembl IDs, use ensembl_gene_id as filter value
go <- getBM(attributes = c("ensembl_gene_id", "go_id"),
            filters = "ensembl_gene_id",
            values = rownames(eset_data_sub),
            mart = ensembl)

# Use the GO.db library to add a column with the GO-term to the dataframe
go$term <- Term(go$go_id)

# Create a named list of character vectors out of the df
s = split(go$ensembl_gene_id, paste(go$go_id,go$term))

# Saves the list as a R environment
go.env <- list2env(s)

# Test
stopifnot(class(go.env) == "environment")

# Pagoda pathway analysis ------------------------------------------------------
cat(sprintf("\n\nPagoda pathway analysis...\n\n"))

fname_pagoda <- paste0("../data/pwpca-eset-dep", append, ".rds")

if (file.exists(fname_pagoda)) {
  pwpca_eset_dep <- readRDS(fname_pagoda)
} else {
  set.seed(1106)
  pwpca_eset_dep <- pagoda.pathway.wPCA(varinfo_dep, go.env, n.components = 1,
                                      n.cores = 8)

  saveRDS(pwpca_eset_dep, file = fname_pagoda)
}

# Pagoda gene clusters ---------------------------------------------------------
cat(sprintf("\n\nPagoda gene clusters...\n\n"))

fname_clusters <- paste0("../data/clpca-eset-dep", append, ".rds")

if (file.exists(fname_clusters)) {
  clpca_eset_dep <- readRDS(fname_clusters)
} else {
  clpca_eset_dep <- pagoda.gene.clusters(varinfo_dep,
                                         trim = 7.1 / ncol(varinfo_dep$mat),
                                         n.clusters = 50,
                                         n.cores = 8,
                                         plot = FALSE)

  saveRDS(clpca_eset_dep, file = fname_clusters)
}
