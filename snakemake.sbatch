#!/bin/bash

# sbatch submission script to run main snakemake process. It then submits
# individual jobs from the compute node.

#SBATCH --job-name=snakemake
#SBATCH --output=snakelog.out
#SBATCH --time=24:00:00
#SBATCH --partition=gilad
#SBATCH --mem=4G
#SBATCH --tasks-per-node=4

source activate scqtl

bash submit-snakemake.sh $*
