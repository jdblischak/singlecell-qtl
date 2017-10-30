#!/usr/bin/env Rscript

# Utility script to combine single cell annotation to a tab-delimited file. Also
# outputs the descriptions of the variables in a separate file.

suppressPackageStartupMessages(library("Biobase"))

fname <- Sys.glob("../data/eset/*.rds")
eset <- Reduce(combine, Map(readRDS, fname))

# Annotation
anno <- pData(eset)
write.table(anno, file = "../data/scqtl-annotation.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Annotation descriptions
anno_descrip <- varMetadata(eset)
anno_descrip <- cbind(rownames(anno_descrip), anno_descrip$labelDescription)
colnames(anno_descrip) <- c("variable", "description")
write.table(anno_descrip, file = "../data/scqtl-annotation-description.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
