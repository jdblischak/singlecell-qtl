# Data

* `eset/` - ExpressionSet object per C1 chip saved as RDS files. Versioned with
  Git LFS.

* `lab-info/` - Experimental variables recorded during sample processing. These
  are added to the phenotype data of the ExpressionSet objects in `eset/`. Use
  `pData()` to access this information from R.

* `md5-core/` - md5 checksums for each FASTQ file provided by the core. Used by
  `code/fastq-download.py`.

* `singlecell-qtl.txt.gz` - A gzipped tab-delimited file with counts for human
  genes. Version with Git LFS. Created by `code/output-exp-mat.R`. Make edits to
  the upstream script, not this file.

* `scqtl-annotation.txt` - A tab-delimited file with the annotation for all
  single cells. Created by `code/output-annotation.R`. Instead of editing
  manually, make edits to the upstream script or to `create-expressionset.R`.

* `scqtl-annotation-description.txt` - A tab-delimited file with the column
  descriptions for `scqtl-annotation.txt`. Created by
  `code/output-annotation.R`. Instead of editing manually, make edits to the
  upstream script or to `create-expressionset.R`.

* `batch[1-5]_qc.txt` - Original recording of experimental variables. These are
  processed by `code/scratch-batch[1-5].R` to create the file in
  `lab-info/`. Use the files in `lab-info/`

## Description of how lab information of single cells is processed

1. PoYuan manually enters information into one Excel file per batch, e.g.
`batch1_qc.xlsx`. Because these are binary files, they are not committed to the
Git repository.

1. John manually exports the Excel data to a tab-delimited text file, e.g.
`batch1_qc.txt` and commits it (use `git add -f` to force add because `data/` is
in `.gitignore`).
    a. After exporting from Excel, fix the line endings by importing into R with
    `readLines` and exporting with `writeLines`.
    b. The script should be executed from `code/`.

1. John manually runs a small script, e.g. `code/scratch-batch1.R`, to generate
one QC file per C1 chip in `data/lab-info/`.
    a. Batch 1 is so different from the other batches (spike-ins and sample
    swaps) that it wasn't worth building one universal script for this step.

1. The Snakemake pipeline automatically detects the files in `data/lab-info/`
and executes the script `code/create-expressionset.R`. This inserts the lab-info
into the ExpressionSet objects in `data/eset/`. The information is then
accessible from the ExpressionSet object using `pData`.
