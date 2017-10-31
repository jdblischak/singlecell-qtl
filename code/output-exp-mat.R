#!/usr/bin/env Rscript

# Utility script to export counts for human genes to gzipped tab-delimited file.
#
# Usage:
#
#   Rscript output-exp-mat.R <directory with ExpressionSet rds files>
#                            <output gzipped tab-delimited file>
#

suppressPackageStartupMessages(library("Biobase"))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
dir_eset <- args[1]
outfile <- args[2]

fname <- Sys.glob(file.path(dir_eset, "*.rds"))
eset <- Reduce(combine, Map(readRDS, fname))

# Only keep human genes
eset_hs <- eset[fData(eset)$source == "H. sapiens", ]

# Add individual to column names
mat_hs <- exprs(eset_hs)
colnames(mat_hs) <- paste(eset_hs$chip_id, eset_hs$experiment, eset_hs$well,
                          sep = ".")

# Make rownames of gene IDs a separate column to make it easier to read with
# data.table::fread.
mat_df <- data.frame(gene = rownames(mat_hs), mat_hs)

gz_outfile <- gzfile(outfile, "w")
write.table(mat_df, file = gz_outfile, quote = FALSE, sep = "\t",
            row.names = FALSE)
close(gz_outfile)
