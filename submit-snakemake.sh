#!/bin/bash

source activate scqtl

snakemake \
    -kp \
    --ri \
    -j 500 \
    --max-jobs-per-second 5 \
    --cluster-config cluster.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --partition=broadwl \
        --job-name={cluster.name} \
	--output={cluster.logfile}" \
    $*
