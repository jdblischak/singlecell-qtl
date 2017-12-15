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
