# Microbial RNAseq
Analysis of RNAseq data from (host-associated) microbial mixtures

### Analysis Workflow

  1. Reads are aligned against the host genome (e.g. human) and all aligning reads are removed.

  2. Non-host reads are aligned against all of the rRNA genes from a set of genomes.

  3. All of the genomes with >= X% coverage per sample, across the rRNA genes, are selected.

  4. Non-host reads are aligned against the set of genomes which have >X% in the ribosomes.

  5. Summary statistics are computed for every genome.


### Input Files

#### Batch file

A group of samples is defined in a comma-delimited "batch file", which includes a single line for each sample.
The name of each sample is given in the column with the header `name`. Samples with single-ended sequencing
or paired-end sequencing which has been interleaved will use the column `fastq` to point to the FASTQ file
(gzip optional) with that data. Paired-end sequencing experiments in which the data is available in two
files (one for the forward read and one for the reverse read) will use the columns `fastq1` and `fastq2`.

#### Reference Genomes

Each reference genome has:

  * A unique name (probably best to avoid any whitespaces)
  * A gzipped FASTA file
  * A gzipped GFF3 file

For example (`test/genome_list.csv`):

```
Escherichia_coli_str_K-12_substr_MG1655,ref_genomes/GCF_000005845.2_ASM584v2_genomic.fna.gz,ref_genomes/GCF_000005845.2_ASM584v2_genomic.gff.gz
Bacteroides_fragilis_YCH46,ref_genomes/GCF_000009925.1_ASM992v1_genomic.fna.gz,ref_genomes/GCF_000009925.1_ASM992v1_genomic.gff.gz
Shigella_sonnei_53G,ref_genomes/GCF_000283715.1_ASM28371v1_genomic.fna.gz,ref_genomes/GCF_000283715.1_ASM28371v1_genomic.gff.gz
```

These reference genomes can be concatenated and indexed using the `build_database.nf` workflow.


### Example database

We have compiled an example reference database and it is currently hosted at 
`s3://fh-ctr-public-reference-data/tool_specific_data/microbial-rnaseq/2019-07-03/`.
An example of running input data against this database is given below.


#### Host Genome

In addition, you must specify a host genome (e.g. human) which will be used to 
subtract data which shouldn't be aligned against bacteria.

This action is also performed with the `build_database.nf` workflow.


### Output files

The output of the tool is specified with two flags:

  * `--output_folder`: Location to place all of the data
  * `--output_prefix`: String that will form the beginning of all output files


### Parameter(s)

At the moment, this tool only takes one tunable parameter: `--min_cov_pct`, which is
the minimum amount of coverage across the ribosomes that a genome must have in order 
to justify a full genome alignment.


### Running the tool

The tool can be run directly from this GitHub repo, e.g.:

```

# You will probably want to change all of these parameters to run on your dataset.

nextflow \
    run \
    fredhutch/microbial-rnaseq \
    --batchfile batchfile.csv \
    --interleaved \
    --host_genome "s3://fh-ctr-public-reference-data/tool_specific_data/microbial-rnaseq/2019-06-10/Homo_sapiens_assembly38.fasta.tar" \
    --database_folder "s3://fh-ctr-public-reference-data/tool_specific_data/microbial-rnaseq/2019-07-03/" \
    --database_prefix 2019-07-03-rnaseq-database \
    --min_cov_pct 10 \
    --output_folder results/ \
    --output_prefix 2019-05-08-test \
    -work-dir work/ \
    -resume

```


### Outputs

This tool creates many CSV files in the `--output_folder`. These include:

  * `<output_prefix>.mapping_summary.csv`: Number of reads input, filtered against host, aligned against genomes.
  * `<output_prefix>.all.csv`: Depth of sequencing for every gene across all samples in all organisms (long format).
  * `<output_prefix>.CDS.csv`: The weighted average depth of sequencing across all CDS genes, broken out by organism and sample.
  * `<output_prefix>.rRNA.csv`: The weighted average depth of sequencing across all rRNA genes, broken out by organism and sample.
  * `<output_prefix>.rRNA_CDS_ratio.csv`: The weighted average depth of sequencing of the rRNA genes / CDS genes, broken out by organism and sample.
  * `<output_prefix>.<organism_name>.csv`: The depth of sequencing for each gene in a given organism, across all samples


### Reference Genomes

You can get a set of reference genomes from NCBI by going to 
https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/ and selecting those genomes with (for example)
an assembly level of "Complete" listed. I made a small script to convert that list into the explicit
format required as input for this pipeline: `convert_ncbi_genome_list.py`.
