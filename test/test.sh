#!/bin/bash

set -e

nextflow \
    -C nextflow.config \
    run \
    ../main.nf \
    --batchfile batchfile.csv \
    --genome_list genome_list.csv \
    --host_genome ref_genomes/AY064377.1.fasta \
    --min_cov_pct 10 \
    --output_folder results/ \
    --output_prefix 2019-05-08-test \
    -with-docker ubuntu:16.04 \
    -process.executor 'local' \
    -work-dir work/ \
    -resume
