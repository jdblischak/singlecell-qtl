#!/usr/bin/env Rscript

# Prepare the metadata to copy-paste in the Excel spreadsheet form for
# submission to GEO.
#
# This file was run after running code/geo.R via code/geo.sh.
#
# Creates 3 files:
# * geo-samples.txt for the section "Samples"
# * geo-processed.txt for the section "Processed Data Files"
# * geo-raw.txt for the section "Raw Files"
#
# A copy of these files are saved both in the directory with the FASTQ files to
# be uploaded to GEO as well as the Git repository to be versioned.

library("stringr")
library("tools")

# Hard-coded variables ---------------------------------------------------------

# dir_geo was created by code/geo.R
dir_geo <- "/project2/gilad/singlecell-qtl/geo"
repo <- "/home/jdblischak/singlecell-qtl"

# Input ------------------------------------------------------------------------

fq <- list.files(path = dir_geo, pattern = "fastq.gz$")
# Read in annotation to get identification
anno <- read.delim(file.path(repo, "data/scqtl-annotation.txt"),
                   stringsAsFactors = FALSE)
# Add back leading zero for experiment ID to get MMDDYYYY
anno$experiment <- sprintf("%08d", anno$experiment)

# Samples ----------------------------------------------------------------------

samples <- list()
samples$name <- sprintf("Sample %04d", seq_along(fq))
samples$title <- str_replace(fq, ".fastq.gz$", "")
samples$source <- rep("LCL-derived iPSC", length(fq))
samples$organism <- rep("Homo sapiens", length(fq))
samples$experiment <- str_extract(fq, "^[0-9]{8}")
samples$well <- str_extract(fq, "[A-H][0-1][0-9]")
samples$individual <- str_extract(fq, "NA[0-9]{5}")
# Add batch from annotation file
stopifnot(anno$experiment == samples$experiment,
          anno$well == samples$well,
          anno$chip_id == samples$individual)
samples$batch <- anno$batch
samples$description <- rep("", length(fq))
samples$molecule <- rep("polyA RNA", length(fq))
# processed data files
samples$processed1 <- rep("eset.rds", length(fq))
samples$processed2 <- rep("scqtl-counts.txt.gz", length(fq))
samples$processed3 <- rep("scqtl-annotation.txt", length(fq))
samples$processed4 <- rep("scqtl-annotation-description.txt", length(fq))
samples$processed5 <- rep("quality-single-cells.txt", length(fq))
samples$processed6 <- rep("genes-pass-filter.txt", length(fq))
samples$raw <- fq
samples_df <- as.data.frame(samples)
head(samples_df)
# Save to GEO directory for upload
samples_fname <- file.path(dir_geo, "geo-samples.txt")
write.table(samples_df, file = samples_fname,
            quote = FALSE, sep = "\t", row.names = FALSE)
# Save to Git repo to version
samples_fname_git <- file.path(repo, "data/geo-samples.txt")
write.table(samples_df, file = samples_fname_git,
            quote = FALSE, sep = "\t", row.names = FALSE)

# Processed data files ---------------------------------------------------------

processed <- list()
processed$name <- c("eset.rds",
                    "scqtl-counts.txt.gz",
                    "scqtl-annotation.txt",
                    "scqtl-annotation-description.txt",
                    "quality-single-cells.txt",
                    "genes-pass-filter.txt")
processed$type <- c("Serialized R object of Bioconductor ExpressionSet class",
                    "gzipped tab-delimited text",
                    rep("tab-delimited text", 4))
stopifnot(file.exists(file.path(repo, "data", processed$name)))
processed$checksum <- md5sum(file.path(repo, "data", processed$name))
processed_df <- as.data.frame(processed)
processed_df
# Save to GEO directory for upload
processed_fname <- file.path(dir_geo, "geo-processed.txt")
write.table(processed_df, file = processed_fname,
            quote = FALSE, sep = "\t", row.names = FALSE)
# Save to Git repo to version
processed_fname_git <- file.path(repo, "data/geo-processed.txt")
write.table(processed_df, file = processed_fname_git,
            quote = FALSE, sep = "\t", row.names = FALSE)

# Raw files --------------------------------------------------------------------

raw <- list()
raw$name <- fq
raw$type <- rep("fastq", length(fq))
# MD5 checksums were computed with code/geo.R
md5_fname <- list.files(path = dir_geo, pattern = ".md5$", full.names = TRUE)
md5_checksums <- Map(readLines, md5_fname)
md5_checksums <- as.character(md5_checksums)
raw$checksum <- md5_checksums
raw$instrument <- rep("HiSeq 2500", length(fq))
raw$read_length <- rep("50", length(fq))
raw$library <- rep("single", length(fq))
raw_df <- as.data.frame(raw)
head(raw_df)
# Save to GEO directory for upload
raw_fname <- file.path(dir_geo, "geo-raw.txt")
write.table(raw_df, file = raw_fname,
            quote = FALSE, sep = "\t", row.names = FALSE)
# Save to Git repo to version
raw_fname_git <- file.path(repo, "data/geo-raw.txt")
write.table(raw_df, file = raw_fname_git,
            quote = FALSE, sep = "\t", row.names = FALSE)
