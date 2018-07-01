#!/usr/bin/env Rscript

# Create exons file for mapping reads to genes with featureCounts.
#
# Usage:
#   Rscript create-exons.R archive org [chr1 chr2 ...] > file.saf
#
# e.g.
#   Rscript create-exons.R feb2014.archive.ensembl.org hs 1 2 3 4 5 > file.saf
#
# where archive is an Ensembl archive URL and org is a 2 letter abbreviation for
# the organism.
#
# Notes:
# + This includes only coding genes, i.e. gene_biotype == "protein_coding" and
#   transcript_biotype == "protein_coding"
# + Output is in Simplified Annotation Format (SAF)
#     + Columns: GeneID, Chr, Start, End, Strand
#     + Coordinates are 1-based, inclusive on both ends
# + Contains duplicate and overlapping exons (featureCounts handles this)

# Parse input arguments
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
archive <- args[1]
stopifnot(grepl("ensembl", archive))
org <- args[2]
stopifnot(org %in% c("ce", "dm", "hs"))
chroms <- args[-1:-2]

# Determine organism dataset
if (org == "ce") {
  dataset <- "celegans_gene_ensembl"
} else if (org == "dm") {
  dataset <- "dmelanogaster_gene_ensembl"
} else if (org == "hs") {
  dataset <- "hsapiens_gene_ensembl"
} else {
  stop(sprintf("Unrecognized organism abbreviation: %s", org))
}

# Download exons from Ensembl database
suppressMessages(library("biomaRt"))
ensembl <- useMart(host = archive,
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = dataset)
# The name of the attribute to obtain the gene name has changed over different
# Ensembl versions (I didn't try to pinpoint the exact switch point). Thus need
# to obtain this dynamically based on the input.
#
# Ensembl 75 (GRCh37): external_gene_id "Associated Gene Name"
#
# Ensembl 88 (GRCh38): external_gene_name "Gene name"
atts <- listAttributes(ensembl, page = "feature_page")
potential_gene_name <- union(match("Associated Gene Name", atts$description),
                             match("Gene name", atts$description))
potential_gene_name <- potential_gene_name[!is.na(potential_gene_name)]
att_gene_name <- atts[potential_gene_name, "name"]
stopifnot(is.character(att_gene_name), length(att_gene_name) == 1)
exons_all <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id",
                                  "chromosome_name", "exon_chrom_start",
                                  "exon_chrom_end", "strand",
                                  att_gene_name,
                                  "gene_biotype", "transcript_biotype"),
                   mart = ensembl)
# Filter by chromosomes and biotype (only include protein_coding)
exons_final <- exons_all[exons_all$chromosome_name %in% chroms &
                         exons_all$gene_biotype == "protein_coding" &
                         exons_all$transcript_biotype == "protein_coding",
                         c("ensembl_gene_id", "chromosome_name", "exon_chrom_start",
                           "exon_chrom_end", "strand", att_gene_name)]
colnames(exons_final) <- c("GeneID", "Chr", "Start", "End", "Strand", "Name")
# Sort by chromosome and position
exons_final <- exons_final[order(exons_final$Chr,
                                 exons_final$Start,
                                 exons_final$End), ]
# Fix chromosome names
exons_final$Chr <- paste0(org, exons_final$Chr)
# Fix strand
exons_final$Strand <- ifelse(exons_final$Strand == 1, "+", "-")

# Save as tab-separated file in Simplified Annotation Format (SAF)
write.table(exons_final, "", quote = FALSE, sep = "\t",
            row.names = FALSE)
