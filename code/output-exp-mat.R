#!/usr/bin/env Rscript

# Utility script to export counts for human genes to gzipped tab-delimited file.

suppressPackageStartupMessages(library("Biobase"))

fname <- Sys.glob("../data/eset/*.rds")
eset <- Reduce(combine, Map(readRDS, fname))

# Only keep cells with 1 cell per well, TRA-1-60 staining for pluripotency, and
# a valid assignment for individual (must be 1 of the 4 individuals on that
# chip)
eset_quality <- eset[, eset$cell_number == 1 & eset$tra1.60 & eset$valid_id]

# Only keep human genes
eset_hs <- eset_quality[fData(eset_quality)$source == "H. sapiens", ]

# Add individual to column names
mat_hs <- exprs(eset_hs)
colnames(mat_hs) <- paste(eset_hs$chip_id, eset_hs$experiment, eset_hs$well,
                          sep = ".")

out <- gzfile("../data/singlecell-qtl.txt.gz", "w")
write.table(mat_hs, file = out, quote = FALSE, sep = "\t",
            row.names = FALSE)
close(out)
