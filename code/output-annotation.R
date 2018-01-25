#!/usr/bin/env Rscript

# Utility script to combine single cell annotation to a tab-delimited file. Also
# outputs the descriptions of the variables in a separate file.
#
# Usage:
#
#   Rscript output-annotation.R <ExpressionSet rds file>
#                               <output annotation file>
#                               <output description file>
#
suppressPackageStartupMessages(library("Biobase"))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
fname_eset <- args[1]
stopifnot(file.exists(fname_eset))
fname_anno <- args[2]
fname_description <- args[3]

eset <- readRDS(fname_eset)

# Annotation

anno <- pData(eset)
write.table(anno, file = fname_anno,
            quote = FALSE, sep = "\t", row.names = FALSE)

# Annotation descriptions
anno_descrip <- varMetadata(eset)
anno_descrip <- cbind(rownames(anno_descrip), anno_descrip$labelDescription)
colnames(anno_descrip) <- c("variable", "description")
write.table(anno_descrip, file = fname_description,
            quote = FALSE, sep = "\t", row.names = FALSE)
