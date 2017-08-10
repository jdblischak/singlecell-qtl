#!/usr/bin/env Rscript

# Create exons file of ERCC spike-ins for mapping reads to genes with featureCounts.
#
# Usage:
#   Rscript create-exons-ercc.R path/to/ERCC92.gtf > file.saf
#
# where ERCC92.gtf is downloaded from Invitrogen:
# http://media.invitrogen.com.edgesuite.net/softwares/ERCC92.gtf
#
# Notes:
# + Output is in Simplified Annotation Format (SAF)
#     + Columns: GeneID, Chr, Start, End, Strand
#     + Coordinates are 1-based, inclusive on both ends

# Parse input arguments
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
ercc <- args[1]
stopifnot(basename(ercc) == "ERCC92.gtf")

# Import ERCC data
ercc_gtf <- read.table(ercc, sep = "\t", stringsAsFactors = FALSE)
# http://www.genome.ucsc.edu/FAQ/FAQformat.html#format3
colnames(ercc_gtf) <- c("seqname", # The name of the sequence. Must be a chromosome or scaffold.
                        "source",  # The program that generated this feature.
                        "feature", # The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
                        "start",   # The starting position of the feature in the sequence. The first base is numbered 1.
                        "end",     # The ending position of the feature (inclusive).
                        "score",   # A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
                        "strand",  # Valid entries include '+', '-', or '.' (for don't know/don't care).
                        "frame",   # If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
                        "group"    # All lines with the same group are linked together into a single item.
)
ercc_saf <- ercc_gtf[, c("seqname", "seqname", "start", "end", "strand",
                         "seqname")]
colnames(ercc_saf) <- c("GeneID", "Chr", "Start", "End", "Strand", "Name")

# Save as tab-separated file in Simplified Annotation Format (SAF)
write.table(ercc_saf, "", quote = FALSE, sep = "\t",
            row.names = FALSE)
