#!/bin/bash

set -e

nextflow \
    -C nextflow.config \
    run \
    ../main.nf \
    --batchfile batchfile.csv \
    --genome_list genome_list.csv \
    --host_genome ref_genomes/AY064377.1.fasta \
    -with-docker ubuntu:16.04 \
    -process.executor 'local' \
    -work-dir work/ \
    -resume
