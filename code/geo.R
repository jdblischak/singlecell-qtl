#!/usr/bin/env Rscript

# In preparation for GEO submission, perform the following 2 steps:
#
# 1. Create a copy of the combined FASTQ file that includes the assinged
# individual: <experiment>-<well>-<individual>.fastq.gz, e.g.
# 02192018-A01-NA18517.fastq.gz
#
# 2. Calculate the MD5 checksum. Save as
# <experiment>-<well>-<individual>.fastq.gz, e.g. 02192018-A01-NA18517.md5
#
#
# This script is called by geo.sh to parallelize across FASTQ files.

library("stringr")
library("tools")

# Hard-coded variables ---------------------------------------------------------

dir_fq_new <- "/project2/gilad/singlecell-qtl/geo"
dir.create(dir_fq_new, showWarnings = FALSE)
repo <- "/home/jdblischak/singlecell-qtl"

# Input ------------------------------------------------------------------------

# Obtain FASTQ file from command-line
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
fastq <- args[1]
# fastq <- "/project2/gilad/singlecell-qtl/fastq-combined/02192018/02192018-A01.fastq.gz"
stopifnot(file.exists(fastq))

# Read in annotation to get identification
anno <- read.delim(file.path(repo, "data/scqtl-annotation.txt"),
                   stringsAsFactors = FALSE)
# Add back leading zero for experiment ID to get MMDDYYYY
anno$experiment <- sprintf("%08d", anno$experiment)

# Obtain samples attributes ----------------------------------------------------

experiment <- str_extract(basename(fastq), "[0-9]{8}")
well <- str_extract(basename(fastq), "[A-Z][0-9]{2}")
ind <- anno$chip_id[anno$experiment == experiment & anno$well == well]

# Copy file --------------------------------------------------------------------

fastq_copy <- paste(experiment, well, ind, sep = "-")
fastq_copy <- paste0(fastq_copy, ".fastq.gz")
fastq_copy <- file.path(dir_fq_new, fastq_copy)

stopifnot(file.copy(fastq, fastq_copy, overwrite = TRUE))

# MD5 checksum -----------------------------------------------------------------

md5 <- md5sum(fastq_copy)

md5_fname <- paste(experiment, well, ind, sep = "-")
md5_fname <- paste0(md5_fname, ".md5")
md5_fname <- file.path(dir_fq_new, md5_fname)

writeLines(md5, con = md5_fname)
