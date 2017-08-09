# Snakefile
#
# This Snakefile runs the single cell analysis starting from the fastq files.
#
# To configure the paths to data files and other settings, edit
# config.yaml.
#
# To configure job submission settings for your cluster, edit
# cluster.json and submit-snakemake.sh.
#
# To run on RCC Midway2 use `bash submit-snakemake.sh`

import glob
import os
from snakemake.utils import R

# Configuration ----------------------------------------------------------------

configfile: "config.yaml"

# Specify Ensembl release for genome sequence and annotation
ensembl_archive = config["ensembl_archive"]
ensembl_rel = config["ensembl_rel"]
ensembl_ftp = "ftp://ftp.ensembl.org/pub/release-" + \
              str(ensembl_rel) + "/fasta/"
ensembl_exons = "exons-ensembl-release-" + str(ensembl_rel) + ".saf"
ensembl_genome_ce = config["ensembl_genome_ce"]
ensembl_genome_dm = config["ensembl_genome_dm"]
ensembl_genome_hs = config["ensembl_genome_hs"]

# Paths to data (must end with forward slash)
dir_proj = config["dir_proj"]
dir_data = dir_proj + "data/"
dir_fq = dir_data + "fastq/"
scratch = config["scratch"]
dir_genome = scratch + "genome-ensembl-release-" + str(ensembl_rel) + "/"
dir_fq_tmp = scratch + "sstmp-fastq/"
dir_bam = scratch + "sstmp-bam/"
dir_bam_dedup = scratch + "sstmp-bam-dedup/"
dir_bam_dedup_stats = scratch + "sstmp-bam-dedup-stats/"
dir_counts = scratch + "sstmp-counts/"
dir_id = dir_data + "id/"

assert os.path.exists(dir_proj), "Project directory exists"
assert os.path.exists(scratch), "Scratch directory exists"

# Directory to send log files. Needs to be created manually since it
# is not a file created by a Snakemake rule.
dir_log = config["dir_log"]
if not os.path.isdir(dir_log):
    os.mkdir(dir_log)

# Names of chromosomes
chr_ce = config["chr_ce"]
chr_dm = config["chr_dm"]
chr_hs = config["chr_hs"]
    
# Input samples ----------------------------------------------------------------

samples = glob_wildcards(dir_fq + "{samples}.fastq.gz").samples

# Targets ----------------------------------------------------------------------

rule all:
    input: expand(dir_bam_dedup + "{sample}-dedup.bam", sample = samples)

rule test_one:
    input: dir_bam_dedup + "YG-PYT-03172017-D08_S428_L004_R1_001-dedup.bam"

rule test_more:
    input: expand(dir_bam_dedup + "{sample}-dedup.bam", sample = samples[:10])

rule download_fasta:
    input: expand(dir_genome + "Caenorhabditis_elegans." + ensembl_genome_ce + \
                  ".dna_sm.chromosome.{chr}.fa.gz", chr = chr_ce),
           expand(dir_genome + "Drosophila_melanogaster." + ensembl_genome_dm + \
                  ".dna_sm.chromosome.{chr}.fa.gz", chr = chr_dm),
           expand(dir_genome + "Homo_sapiens." + ensembl_genome_hs + \
                  ".dna_sm.chromosome.{chr}.fa.gz", chr = chr_hs),
           dir_genome + "ercc.fa"

# Prepare genome annotation ----------------------------------------------------

rule download_genome_ce:
    output: dir_genome + "Caenorhabditis_elegans." + ensembl_genome_ce + \
            ".dna_sm.chromosome.{chr}.fa.gz"
    params: chr = "{chr}", build = ensembl_genome_ce,
            ftp = ensembl_ftp + "caenorhabditis_elegans/dna/"
    shell: "wget -O {output} {params.ftp}Caenorhabditis_elegans.{params.build}.dna_sm.chromosome.{params.chr}.fa.gz"

rule download_genome_dm:
    output: dir_genome + "Drosophila_melanogaster." + ensembl_genome_dm + \
            ".dna_sm.chromosome.{chr}.fa.gz"
    params: chr = "{chr}", build = ensembl_genome_dm,
            ftp = ensembl_ftp + "drosophila_melanogaster/dna/"
    shell: "wget -O {output} {params.ftp}Drosophila_melanogaster.{params.build}.dna_sm.chromosome.{params.chr}.fa.gz"

rule download_genome_hs:
    output: dir_genome + "Homo_sapiens." + ensembl_genome_hs + \
            ".dna_sm.chromosome.{chr}.fa.gz"
    params: chr = "{chr}", build = ensembl_genome_hs,
            ftp = ensembl_ftp + "homo_sapiens/dna/"
    shell: "wget -O {output} {params.ftp}Homo_sapiens.{params.build}.dna_sm.chromosome.{params.chr}.fa.gz"

rule download_ercc:
    output: dir_genome + "ercc.fa"
    shell: "wget -O {output} http://tools.invitrogen.com/downloads/ERCC92.fa"

rule unzip_chromosome_fasta_ce:
    input: expand(dir_genome + "Caenorhabditis_elegans." + ensembl_genome_ce + \
                  ".dna_sm.chromosome.{chr}.fa.gz", chr = chr_ce)
    output: temp(dir_genome + "ce.fa")
    shell: "zcat {input} | sed 's/>/>ce/' > {output}"

rule unzip_chromosome_fasta_dm:
    input: expand(dir_genome + "Drosophila_melanogaster." + ensembl_genome_dm + \
                  ".dna_sm.chromosome.{chr}.fa.gz", chr = chr_dm)
    output: temp(dir_genome + "dm.fa")
    shell: "zcat {input} | sed 's/>/>dm/' > {output}"

rule unzip_chromosome_fasta_hs:
    input: expand(dir_genome + "Homo_sapiens." + ensembl_genome_hs + \
                  ".dna_sm.chromosome.{chr}.fa.gz", chr = chr_hs)
    output: temp(dir_genome + "hs.fa")
    shell: "zcat {input} | sed 's/>/>hs/' > {output}"

# Quanitify expression with Subjunc/featureCounts ------------------------------

rule subread_index:
    input: dir_genome + "ce.fa",
           dir_genome + "dm.fa",
           dir_genome + "hs.fa",    
           dir_genome + "ercc.fa"
    output: dir_genome + "genome.reads"
    params: prefix = dir_genome + "genome"
    shell: "subread-buildindex -o {params.prefix} {input}"

rule extract_umi:
    input: dir_fq + "{sample}.fastq.gz"
    output: temp(dir_fq_tmp + "{sample}.fastq.gz")
    shell: "umi_tools extract --bc-pattern=NNNNNNNNN -I {input} -S {output}"

rule subjunc:
    input: read = dir_fq_tmp + "{sample}.fastq.gz",
           index = dir_genome + "genome.reads"
    output: temp(dir_bam + "{sample}.bam")
    params: prefix = dir_genome + "genome"
    threads: 8
    shell: "subjunc -i {params.prefix} -r {input.read} -T {threads} > {output}"

rule sort_bam:
    input: dir_bam + "{sample}.bam"
    output: dir_bam + "{sample}-sort.bam"
    shell: "samtools sort -o {output} {input}"

rule index_bam:
    input: dir_bam + "{sample}-sort.bam"
    output: dir_bam + "{sample}-sort.bam.bai"
    shell: "samtools index {input}"

rule dedup_umi:
    input: bam = dir_bam + "{sample}-sort.bam",
           index = dir_bam + "{sample}-sort.bam.bai"
    output: bam = dir_bam_dedup + "{sample}-dedup.bam",
            edit_distance = dir_bam_dedup_stats + "{sample}_edit_distance.tsv",
            per_umi_per_position = dir_bam_dedup_stats + "{sample}_per_umi_per_position.tsv",
            per_umi = dir_bam_dedup_stats + "{sample}_per_umi.tsv"
    params: stats = dir_bam_dedup_stats + "{sample}"
    shell: "umi_tools dedup -I {input.bam} --output-stats={params.stats} -S {output.bam}"
