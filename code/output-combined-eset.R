#!/usr/bin/env Rscript

# Utility script to combine the ExpressionSet RDS files for each chip into a
# combined ExpressionSet RDS object for quick import.
#
# Usage:
#
# Rscript output-combined-eset.R <directory with ExpressionSet rds files> \
#                                <output combined RDS file>
#

suppressPackageStartupMessages(library("Biobase"))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
dir_eset <- args[1]
stopifnot(dir.exists(dir_eset))
outfile <- args[2]

fname <- Sys.glob(file.path(dir_eset, "*.rds"))
eset <- Reduce(combine, Map(readRDS, fname))

saveRDS(eset, file = outfile)
