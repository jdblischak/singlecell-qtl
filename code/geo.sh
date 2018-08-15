#!/bin/bash

# Submit geo.R to prepare FASTQ files for GEO submission

dir="/project2/gilad/singlecell-qtl/fastq-combined"
repo="/home/jdblischak/singlecell-qtl"

for fq in $dir/083*/*fastq.gz
do
 name=`basename $fq`
 echo $name
 sbatch -p broadwl --job-name=$name --mem=1g --wrap="Rscript $repo/code/geo.R $fq"
 sleep 1s
done
