# Snakefile
#
# This Snakefile runs the single cell analysis starting from the FASTQ files.
#
# To configure the paths to data files and other settings, edit
# config.yaml.
#
# To configure job submission settings for your cluster, edit
# cluster.json and submit-snakemake.sh.
#
# To run on RCC Midway2, follow these steps:
#
# 1. Start an interactive session. The example command below starts an
#    interactive session on a Midway2 compute node with 4 GB of RAM, 4 CPUs for
#    multithreading, and a time limit of 24 hours.
#
#    sinteractive --mem=4G --tasks-per-node=4 --partition=broadwl --time=24:00:00
#
# 2. Run the following command from the root of the project directory. nohup
#    prevents the job from dying if you lose connection with the session. The
#    script submit-snakemake.sh activates the conda environment and then runs
#    snakemake with reasonable defaults.
#
#    nohup bash submit-snakemake.sh &
#
# 3. Monitor progress using any of the following:
#
#    tail nohup.out
#    grep % nohup.out | tail
#    ls log/ | wc -l
#

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
dir_fq = config["dir_fq"]
dir_fq_combin = dir_data + "fastq-combined/"
scratch = config["scratch"]
dir_genome = scratch + "genome-ensembl-release-" + str(ensembl_rel) + "/"
dir_fq_filter = scratch + "scqtl-fastq-filter/"
dir_fq_extract = scratch + "scqtl-fastq-extract/"
dir_bam = scratch + "scqtl-bam/"
dir_bam_dedup = scratch + "scqtl-bam-dedup/"
dir_bam_dedup_stats = scratch + "scqtl-bam-dedup-stats/"
dir_bam_verify = scratch + "scqtl-bam-verify/"
dir_counts = scratch + "scqtl-counts/"
dir_totals = scratch + "scqtl-totals/"
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

chips = config["chips"]
rows = config["rows"]
cols = config["cols"]

# Constrain wildcards. Necessary to resolve some of the complex rules.
# https://docs.python.org/3/howto/regex.html
wildcard_constraints: chip = "[0-9]{8,8}", row = "[A-H]", col = "[0-1][0-9]"

# Targets ----------------------------------------------------------------------

rule all:
    input: expand(dir_data + "eset/{chip}.rds", chip = chips)

rule chip_03232017:
    input: dir_data + "eset/03232017.rds"

rule chip_04202017:
    input: dir_data + "eset/04202017.rds"

# Functions --------------------------------------------------------------------

# Find all fastq.gz files for a given sample.
# Inspired by this post on the Snakemake Google Group:
# https://groups.google.com/forum/#!searchin/snakemake/multiple$20input$20files%7Csort:relevance/snakemake/bpTnr7FgDuQ/ybacyom6BQAJ
def merge_fastq(wc):
    pattern = dir_fq + "YG-PYT-{chip}-{row}{col}_S{{s}}_L{{lane}}_R1_001.fastq.gz"
    unknowns = glob_wildcards(pattern.format(chip = wc.chip, row = wc.row,
                                             col = wc.col))
    files = expand(pattern.format(chip = wc.chip, row = wc.row, col = wc.col),
                   zip, s = unknowns.s, lane = unknowns.lane)
    return files

# Prepare genome annotation ----------------------------------------------------

localrules: download_ercc, download_ercc_gtf, gather_exons

rule target_exons:
    input: dir_genome + ensembl_exons

rule target_fasta:
    input: expand(dir_genome + "Caenorhabditis_elegans." + ensembl_genome_ce + \
                  ".dna_sm.chromosome.{chr}.fa.gz", chr = chr_ce),
           expand(dir_genome + "Drosophila_melanogaster." + ensembl_genome_dm + \
                  ".dna_sm.chromosome.{chr}.fa.gz", chr = chr_dm),
           expand(dir_genome + "Homo_sapiens." + ensembl_genome_hs + \
                  ".dna_sm.chromosome.{chr}.fa.gz", chr = chr_hs),
           dir_genome + "ercc.fa"

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

rule create_exons:
    output: dir_genome + "{organism}.saf"
    params: archive = ensembl_archive, organism = "{organism}",
            # Hack to dynamically get the list of chromosomes for each organism
            # https://stackoverflow.com/a/45585380/2483477
            chroms = lambda wildcards: globals()["chr_" + "{organism}".format(**wildcards)]
    shell: "Rscript code/create-exons.R {params.archive} {params.organism} \
            {params.chroms} > {output}"

rule download_ercc_gtf:
    output: dir_genome + "ERCC92.gtf"
    shell: "wget -O {output} http://media.invitrogen.com.edgesuite.net/softwares/ERCC92.gtf"

rule create_exons_ercc:
    input: dir_genome + "ERCC92.gtf"
    output: dir_genome + "ercc.saf"
    shell: "Rscript code/create-exons-ercc.R {input} > {output}"

rule gather_exons:
    input: expand(dir_genome + "{organism}.saf", \
                  organism = ["ce", "dm", "ercc", "hs"])
    output: dir_genome + ensembl_exons
    shell: "cat {input[0]} | grep GeneID > {output}; \
           cat {input} | grep -v GeneID >> {output}"

# Quantify expression with Subjunc/featureCounts -------------------------------

localrules: index_bam, index_bam_dedup

rule target_counts:
    input: bam = expand(dir_counts + "{chip}/{chip}-{row}{col}.txt", \
                        chip = chips, row = rows, col = cols)

rule target_bam:
    input: bam = expand(dir_bam + "{chip}/{chip}-{row}{col}-sort.bam", \
                        chip = chips, row = rows, col = cols),
           index = expand(dir_bam + "{chip}/{chip}-{row}{col}-sort.bam.bai", \
                          chip = chips, row = rows, col = cols)

rule target_fastq:
    input: expand(dir_fq_combin + "{chip}/{chip}-{row}{col}.fastq.gz", \
                  chip = chips, row = rows, col = cols)

rule subread_index:
    input: dir_genome + "ce.fa",
           dir_genome + "dm.fa",
           dir_genome + "hs.fa",    
           dir_genome + "ercc.fa"
    output: dir_genome + "genome.reads"
    params: prefix = dir_genome + "genome"
    shell: "subread-buildindex -o {params.prefix} {input}"

rule combine_fastq:
    input: merge_fastq
    output: dir_fq_combin + "{chip}/{chip}-{row}{col}.fastq.gz"
    shell: "zcat {input} | gzip -c > {output}"

rule filter_umi:
    input: dir_fq_combin + "{chip}/{chip}-{row}{col}.fastq.gz"
    output: temp(dir_fq_filter + "{chip}/{chip}-{row}{col}.fastq.gz")
    run:
        from Bio import SeqIO
        import gzip

        path_in = input[0]
        path_out = output[0]
        handle_in = gzip.open(path_in, "rt")
        handle_out = gzip.open(path_out, "wt")

        fq = SeqIO.parse(handle_in, "fastq")
        for read in fq:
            if read.seq[6:9] == "GGG":
                # Only export reads that have a G in positions 7, 8,
                # and 9
                handle_out.write(read.format("fastq"))
                # Note: I didn't use SeqIO.write here b/c I couldn't
                # find a way to write directly to a gzipped file.
                # https://bioinformatics.stackexchange.com/q/892/1302

        handle_in.close()
        handle_out.close()

rule extract_umi:
    input: dir_fq_filter + "{chip}/{chip}-{row}{col}.fastq.gz"
    output: temp(dir_fq_extract + "{chip}/{chip}-{row}{col}.fastq.gz")
    shell: "umi_tools extract --bc-pattern=NNNNNNNNN -I {input} -S {output}"

rule subjunc:
    input: read = dir_fq_extract + "{chip}/{chip}-{row}{col}.fastq.gz",
           index = dir_genome + "genome.reads"
    output: temp(dir_bam + "{chip}/{chip}-{row}{col}.bam")
    params: prefix = dir_genome + "genome"
    threads: 8
    shell: "subjunc -i {params.prefix} -r {input.read} -T {threads} > {output}"

rule sort_bam:
    input: dir_bam + "{chip}/{chip}-{row}{col}.bam"
    output: dir_bam + "{chip}/{chip}-{row}{col}-sort.bam"
    shell: "samtools sort -o {output} {input}"

rule index_bam:
    input: dir_bam + "{chip}/{chip}-{row}{col}-sort.bam"
    output: dir_bam + "{chip}/{chip}-{row}{col}-sort.bam.bai"
    shell: "samtools index {input}"

rule dedup_umi:
    input: bam = dir_bam + "{chip}/{chip}-{row}{col}-sort.bam",
           index = dir_bam + "{chip}/{chip}-{row}{col}-sort.bam.bai"
    output: bam = temp(dir_bam_dedup + "{chip}/{chip}-{row}{col}.bam"),
            edit_distance = dir_bam_dedup_stats + "{chip}/{chip}-{row}{col}_edit_distance.tsv",
            per_umi_per_position = dir_bam_dedup_stats + "{chip}/{chip}-{row}{col}_per_umi_per_position.tsv",
            per_umi = dir_bam_dedup_stats + "{chip}/{chip}-{row}{col}_per_umi.tsv"
    params: stats = dir_bam_dedup_stats + "{chip}/{chip}-{row}{col}"
    shell: "umi_tools dedup -I {input.bam} --output-stats={params.stats} -S {output.bam}"

# These are copied and modified from the above rules for sorting and indexing a
# BAM file because they the deduplicated BAM files are saved in a separate
# directory. I could put them all in the same directory and change the suffix of
# the filename to distinguish them, but that would result in many files in one
# directory. I had tried combining the rules by have the directory as a
# wildcard. This worked great when run sequentially, but my batch submission
# pipeline inserts the wildcards into the job name and log files so that they
# are interpretable. Inserting a path with forward slashes into the log files
# caused them to run for a second and then instantly die without producing a log
# file.
rule sort_bam_dedup:
    input: dir_bam_dedup + "{chip}/{chip}-{row}{col}.bam"
    output: dir_bam_dedup + "{chip}/{chip}-{row}{col}-sort.bam"
    shell: "samtools sort -o {output} {input}"

rule index_bam_dedup:
    input: dir_bam_dedup + "{chip}/{chip}-{row}{col}-sort.bam"
    output: dir_bam_dedup + "{chip}/{chip}-{row}{col}-sort.bam.bai"
    shell: "samtools index {input}"

rule feauturecounts:
    input: bam = dir_bam + "{chip}/{chip}-{row}{col}-sort.bam",
           dedup = dir_bam_dedup + "{chip}/{chip}-{row}{col}-sort.bam",
           exons = dir_genome + ensembl_exons
    output: dir_counts + "{chip}/{chip}-{row}{col}.txt"
    threads: 8
    shell: "featureCounts -a {input.exons} -F SAF -s 1 --read2pos 5 \
            -T {threads} -o {output} {input.bam} {input.dedup}"

rule gather_counts:
    input: expand(dir_counts + "{{chip}}/{{chip}}-{row}{col}.txt", \
                  row = rows, col = cols)
    output: reads = dir_data + "reads/{chip}.txt.gz",
            molecules = dir_data + "molecules/{chip}.txt.gz"
    run:
        import gzip
        import os

        input.sort()

        reads = gzip.open(output.reads, "wt")
        molecules = gzip.open(output.molecules, "wt")

        # Obtain the gene IDs from the first file
        genes = []
        f1 = open(input[0], "r")
        for line in f1:
            if line[0] == "#" or line[:6] == "Geneid":
                continue
            g = line.strip().split("\t")[0]
            genes.append(g)
        f1.close()

        # Write header
        header = "\t".join(["sample",
                            "experiment",
                            "well"] + genes) + "\n"
        reads.write(header)
        molecules.write(header)

        # Obtain and write gene counts for each sample
        for f in input:
            sample = os.path.basename(f).rstrip(".txt")
            experiment, well = sample.split("-")
            reads.write("\t".join([sample, experiment, well]) + "\t")
            molecules.write("\t".join([sample, experiment, well]) + "\t")
            with open(f, "r") as handle:
                n_reads = [""] * len(genes)
                n_molecules = [""] * len(genes)
                i = 0
                for line in handle:
                    if line[0] == "#" or line[:6] == "Geneid":
                        continue
                    cols = line.strip().split("\t")
                    assert int(cols[6]) >= int(cols[7]), \
                        "Reads greater than or equal to molecules"
                    n_reads[i] = cols[6]
                    n_molecules[i] = cols[7]
                    i += 1
            reads.write("\t".join(n_reads) + "\n")
            molecules.write("\t".join(n_molecules) + "\n")

        reads.close()
        molecules.close()

rule expressionset:
    input: lab = dir_data + "lab-info/{chip}.txt",
           totals = dir_data + "totals/{chip}.txt",
           molecules = dir_data + "molecules/{chip}.txt.gz",
           verify = dir_data + "verify/{chip}.txt",
           saf = dir_genome + ensembl_exons
    output: dir_data + "eset/{chip}.rds"
    shell: "Rscript code/create-expressionset.R {input.molecules} \
                                                {input.lab} \
                                                {input.totals} \
                                                {input.verify} \
                                                {input.saf} \
                                                {output}"

# Calculate total counts -------------------------------------------------------

localrules: gather_totals

rule target_totals:
    input: expand(dir_data + "totals/{chip}.txt", \
                  chip = chips, row = rows, col = cols)

rule count_totals:
    input: fastq = dir_fq_combin + "{chip}/{chip}-{row}{col}.fastq.gz",
           bam = dir_bam + "{chip}/{chip}-{row}{col}-sort.bam",
           bam_index = dir_bam + "{chip}/{chip}-{row}{col}-sort.bam.bai",
           dedup = dir_bam_dedup + "{chip}/{chip}-{row}{col}-sort.bam",
           dedup_index = dir_bam_dedup + "{chip}/{chip}-{row}{col}-sort.bam.bai"
    output: dir_totals + "{chip}/{chip}-{row}{col}.txt"
    run:
        # Count the number of raw reads
        import gzip
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc68
        from Bio.SeqIO.QualityIO import FastqGeneralIterator
        raw = 0
        with gzip.open(input.fastq, "rt") as handle:
            for title, seq, qual in FastqGeneralIterator(handle):
                raw += 1

        # Parse the BAM file to obtain:
        #  * The number of reads with a valid UMI
        #  * The number of mapped (and unmapped) reads
        #  * The number of reads mapped to ce, dm, ercc, hs
        # https://pysam.readthedocs.io/en/stable/api.html#pysam.AlignedSegment
        import pysam
        umi = 0
        unmapped = 0
        mapped = 0
        ce = 0
        dm = 0
        ercc = 0
        hs = 0
        bam = pysam.AlignmentFile(input.bam, "rb")
        for read in bam:
            umi += 1
            if read.is_unmapped:
                unmapped += 1
            else:
                mapped += 1
                ref = read.reference_name
                if ref[:2] == "ce":
                    ce += 1
                elif ref[:2] == "dm":
                    dm += 1
                elif ref[:2] == "hs":
                    hs += 1
                else:
                    ercc += 1
        bam.close()

        # Parse the deduplicated BAM file to obtain:
        #  * The number of molecules
        #  * The number of molecules mapped to ce, dm, ercc, hs
        mol = 0
        mol_ce = 0
        mol_dm = 0
        mol_ercc = 0
        mol_hs = 0
        dedup = pysam.AlignmentFile(input.dedup, "rb")
        for read in dedup:
            mol += 1
            ref = read.reference_name
            if ref[:2] == "ce":
                mol_ce += 1
            elif ref[:2] == "dm":
                mol_dm += 1
            elif ref[:2] == "hs":
                mol_hs += 1
            else:
                mol_ercc += 1
        dedup.close()

        # Consistency checks
        assert umi <= raw, \
            "Reads with a UMI less than or equal to raw reads"
        assert mapped + unmapped == umi, \
            "Mapped and unmapped reads sum to reads with a UMI"
        assert ce + dm + hs + ercc == mapped, \
            "Reads mapped to specific genomes sum to mapped reads"
        assert mol < mapped, \
            "Molecules less than reads."
        assert mol_ce + mol_dm + mol_ercc + mol_hs == mol, \
            "Molecules mapped to specific genomes sum to molecules"

        # Export total counts
        outfile = open(output[0], "w")
        outfile.write("\t".join([str(raw),
                                 str(umi),
                                 str(mapped),
                                 str(unmapped),
                                 str(ce),
                                 str(dm),
                                 str(ercc),
                                 str(hs),
                                 str(mol),
                                 str(mol_ce),
                                 str(mol_dm),
                                 str(mol_ercc),
                                 str(mol_hs)]
                      ) + "\n")

rule gather_totals:
    input: expand(dir_totals + "{{chip}}/{{chip}}-{row}{col}.txt", \
                  row = rows, col = cols)
    output: dir_data + "totals/{chip}.txt"
    run:
        import os

        input.sort()

        outfile = open(output[0], "w")
        header = "\t".join(["sample",
                            "experiment",
                            "well",
                            "raw",
                            "umi",
                            "mapped",
                            "unmapped",
                            "reads_ce",
                            "reads_dm",
                            "reads_ercc",
                            "reads_hs",
                            "molecules",
                            "mol_ce",
                            "mol_dm",
                            "mol_ercc",
                            "mol_hs"]) + "\n"
        outfile.write(header)

        for f in input:
            sample = os.path.basename(f).rstrip(".txt")
            experiment, well = sample.split("-")
            outfile.write("\t".join([sample, experiment, well]) + "\t")
            with open(f, "r") as handle:
                outfile.write(handle.read())
        outfile.close()

# Identify individuals with verifyBamID ----------------------------------------

localrules: index_bam_verify, parse_verify, combine_verify

rule target_verify:
    input: expand(dir_data + "verify/{chip}.txt", chip = chips)

# The VCF file doesn't contain mitochondrial SNPs, thus it is easy to
# convert from UCSC chromsome names to Ensembl by just removing the
# "chr".
rule prepare_genos:
    input: dir_data + "snps.hg19.exons.vcf.gz"
    output: dir_data + "snps.grch37.exons.vcf.gz"
    shell: "zcat {input} | sed 's/chr//g' | gzip -c > {output}"

# Prepare BAM files. verifyBamID only accepts chromsome names such as
# "1" or "chr1" so that it can filter non-autosomal chromsomes. Mine
# are "hs1" to distinguish from the other species.
rule prepare_bam:
    input: bam = dir_bam_dedup + "{chip}/{chip}-{row}{col}-sort.bam",
           index = dir_bam_dedup + "{chip}/{chip}-{row}{col}-sort.bam.bai"
    output: bam = temp(dir_bam_verify + "{chip}/{chip}-{row}{col}-sort.bam")
    shell: "samtools view -H {input.bam} | \
            sed -e 's/SN:hs/SN:/g' | \
            samtools reheader - {input.bam} > {output.bam}"

# Copied and modified from earlier rules. See note above for explanation.
rule sort_bam_verify:
    input: dir_bam_verify + "{chip}/{chip}-{row}{col}.bam"
    output: dir_bam_verify + "{chip}/{chip}-{row}{col}-sort.bam"
    shell: "samtools sort -o {output} {input}"

rule index_bam_verify:
    input: dir_bam_verify + "{chip}/{chip}-{row}{col}-sort.bam"
    output: dir_bam_verify + "{chip}/{chip}-{row}{col}-sort.bam.bai"
    shell: "samtools index {input}"

# Run verifyBamID to obtain the best individual match for the BAM file
rule verify_bam:
    input: vcf =  dir_data + "snps.grch37.exons.vcf.gz",
           bam = dir_bam_verify + "{chip}/{chip}-{row}{col}-sort.bam",
           index = dir_bam_verify + "{chip}/{chip}-{row}{col}-sort.bam.bai"
    output: bestSM = dir_id + "{chip}/{chip}-{row}{col}.bestSM",
            depthSM = dir_id + "{chip}/{chip}-{row}{col}.depthSM",
            selfSM = temp(dir_id + "{chip}/{chip}-{row}{col}.selfSM"),
            log = dir_id + "{chip}/{chip}-{row}{col}.log"
    params: prefix = dir_id + "{chip}/{chip}-{row}{col}"
    shell: "verifyBamID --vcf {input.vcf} --bam {input.bam} --best --ignoreRG --out {params.prefix}"

# Parse the various verifyBamID output files into one tab-separated file. Each
# sample has one file with a header row and a data row.
rule parse_verify:
    input: bestSM = dir_id + "{chip}/{chip}-{row}{col}.bestSM",
           depthSM = dir_id + "{chip}/{chip}-{row}{col}.depthSM"
    output: dir_id + "{chip}/{chip}-{row}{col}-results.txt"
    params: id = "{chip}-{row}{col}"
    run:
        bestSM = open(input.bestSM, "rt")
        depthSM = open(input.depthSM, "rt")
        results = open(output[0], "w")

        for line in bestSM:
            # Confirm the header columns
            if line[0] == "#":
                cols = line.strip().split("\t")
                assert cols[0] == "#SEQ_ID" and cols[2] == "CHIP_ID" and \
                       cols[3] == "#SNPS" and cols[4] == "#READS" and \
                       cols[5] == "AVG_DP" and cols[6] == "FREEMIX" and \
                       cols[11] == "CHIPMIX", "bestSM columns are as expected"
            else:
                cols = line.strip().split("\t")
                chip_id = cols[2]
                snps = cols[3]
                reads = cols[4]
                avg_dp = cols[5]
                freemix = cols[6]
                chipmix = cols[11]

        # Report the number of SNPs that had more than a minimum read depth
        depth_min = 1
        snps_w_min = 0
        for line in depthSM:
            # Confirm the header columns
            if line[0] == "#":
                cols = line.strip().split("\t")
                assert cols[1] == "DEPTH" and cols[2] == "#SNPs", \
                       "depthSM columns are as expected"
            else:
                cols = line.strip().split("\t")
                depth = int(cols[1])
                n_snps = int(cols[2])
                if depth >= depth_min:
                    snps_w_min = snps_w_min + n_snps

        out_header = ["sample", "chip_id", "chipmix", "freemix",
                      "snps", "reads", "avg_dp",
                      "min_dp", "snps_w_min"]
        out_cols = [params.id, chip_id, chipmix, freemix,
                    snps, reads, avg_dp,
                    str(depth_min), str(snps_w_min)]
        results.write("\t".join(out_header) + "\n")
        results.write("\t".join(out_cols) + "\n")

        bestSM.close()
        depthSM.close()
        results.close()

# Combine all the samples for a given chip into one file.
rule combine_verify:
    input: expand(dir_id + "{{chip}}/{{chip}}-{row}{col}-results.txt", \
                  row = rows, col = cols)
    output: dir_data + "verify/{chip}.txt"
    shell:
        "head -n 1 {input[0]} > {output};"
        "cat {input} | grep -v \"id\" | sort -k1 >> {output}"
