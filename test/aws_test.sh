#!/bin/bash

set -e

echo """batchfile.bam.csv bam
batchfile.csv single
batchfile.csv interleaved
batchfile.paired.csv paired""" | while read batchfile mode; do

    echo "Processing $batchfile in $mode mode"

    nextflow \
        -C ~/nextflow.config \
        run \
        ../main.nf \
        --batchfile $batchfile \
        --$mode \
        --host_genome "s3://fh-ctr-public-reference-data/tool_specific_data/microbial-rnaseq/2019-06-10/Homo_sapiens_assembly38.fasta.tar" \
        --database_folder "s3://fh-ctr-public-reference-data/tool_specific_data/microbial-rnaseq/2019-06-10/" \
        --database_prefix 2019-05-12-rnaseq-database \
        --output_folder test_output \
        --output_prefix mode_$mode \
        --min_cov_pct 5 \
        -resume
done