#!/usr/bin/env nextflow

// Default values for boolean flags
params.single = false
params.paired = false
params.interleaved = false
params.bam = false
params.help = false

def helpMessage() {
    log.info"""
    Usage:

    nextflow run fredhutch/microbial-rnaseq <ARGUMENTS>
    
    Arguments:
      --batchfile                   CSV file listing samples to analyze (see below)
      --database_folder             Path to folder containing database created with build_database.nf
      --database_prefix             Name of database (made with build_database.nf)
      --host_genome                 Path to TAR file containing host genome (made with build_database.nf)
      --output_folder               Folder to place outputs
      --output_prefix               Name for output files

    Options:
      --min_cov_pct                 Minimum coverage of rRNA sequences used for filtering per-genome (default 90)
      --single                      Input data is a single FASTQ containing unpaired reads (default)
      --paired                      Input data is paired-end FASTQ in two files (otherwise treat as single-ended)
      --interleaved                 Input data is paired-end FASTQ in a single file (otherwise treat as single-ended)
      --bam                         Input data is in BAM format. Compatible with --interleaved for paired-end data.

    Batchfile:
      The batchfile is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `name`. 
      If data is --single or --interleaved, reads are specified by `fastq`.
      If data is --paired, reads are specified by two columns, `fastq1` and `fastq2`.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// --database_prefix is the name for the database files
params.database_prefix = "microbial_genomes"

// --database_folder is the location of the database
params.database_folder = "./"

// The database contains a TAR with the indexed ribosome sequences
params.ribosome_tar = "${params.database_folder}${params.database_prefix}.ribosomes.tar"
ribosome_tar = file("${params.ribosome_tar}")
params.ribosome_tsv = "${params.database_folder}${params.database_prefix}.ribosomes.tsv"
ribosome_tsv = file("${params.ribosome_tsv}")

// The database contains a FASTA and a TSV linking headers to organisms
params.genome_fasta = "${params.database_folder}${params.database_prefix}.fasta.gz"
genome_fasta = file("${params.genome_fasta}")
params.genome_tsv = "${params.database_folder}${params.database_prefix}.tsv.gz"
genome_tsv = file("${params.genome_tsv}")
params.genome_table = "${params.database_folder}${params.database_prefix}.tsv.gz"
genome_table = file("${params.genome_table}")
params.all_gff = "${params.database_folder}${params.database_prefix}.gff.gz"
all_gff = file("${params.all_gff}")

// Minimum coverage of a ribosomal subunit needed for full genome alignment
params.min_cov_pct = 90
params.min_qual = 50

// --host_genome is a TAR containing an indexed FASTA
params.host_genome = "hg38.fasta.tar"
host_genome_tar = file(params.host_genome)

// --output_folder is the folder in which to place the results
params.output_folder = "./"

// --output_prefix is the name to prepend to all output files
params.output_prefix = "microbial-rnaseq"

// Logic to handle different types of input data
if ( params.paired ){
  assert !params.single && !params.interleaved: "--paired cannot be specified with either --single or --interleaved"

  Channel.from(file(params.batchfile))
        .splitCsv(header: true, sep: ",")
        .map { sample ->
        [sample.name, file(sample.fastq1), file(sample.fastq2)]}
        .set{ interleave_ch }

  process interleave {
    container "ubuntu:16.04"
    cpus 1
    memory "2 GB"
    errorStrategy "retry"

    input:
    set sample_name, file(fastq1), file(fastq2) from interleave_ch

    output:
    set sample_name, file("${fastq1}.interleaved.fastq.gz") into concatenate_ch

    """
    set -e

    # Some basic checks that the files exist and the line numbers match
    [[ -s "${fastq1}" ]]
    [[ -s "${fastq2}" ]]
    (( \$(gunzip -c ${fastq1} | wc -l) == \$(gunzip -c ${fastq2} | wc -l) ))

    # Now interleave the files
    paste <(gunzip -c ${fastq1}) <(gunzip -c ${fastq2}) | paste - - - - | awk -v OFS="\\n" -v FS="\\t" '{print(\$1,\$3,\$5,\$7,\$2,\$4,\$6,\$8)}' | gzip -c > "${fastq1}.interleaved.fastq.gz"
    """
      
  }

}
else if (params.bam){
  assert !params.paired: "--bam cannot be specified with --paired (use --single [default] or --interleaved)"
  if (params.single && paired.interleaved){
    assert false: "Cannot specify both --single and --interleaved"
  }

  Channel.from(file(params.batchfile))
        .splitCsv(header: true, sep: ",")
        .map { sample ->
        [sample.name, file(sample.bam)]}
        .set{ bam_ch }

  // Convert the input BAM file into FASTQ format
  process bamToFastq {
    container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
    cpus 1
    memory "2 GB"

    input:
    set sample_name, file(bam) from bam_ch
    
    output:
    set sample_name, file("${bam}.fastq.gz") into concatenate_ch

    afterScript "rm *"

    """
#!/bin/bash

set -e

samtools fastq "${bam}" | gzip -c > "${bam}.fastq.gz"
      """

  }
}
else {
  if ( params.interleaved ){
    assert !params.single && !params.paired: "--interleaved cannot be specified with either --single or --paired"
  }
  if ( params.single ){
    assert !params.interleaved && !params.paired: "--single cannot be specified with either --interleaved or --paired"
  }

  // For either `interleaved` or `single` (including no flags), just put the FASTQ into the same analysis queue
  Channel.from(file(params.batchfile))
        .splitCsv(header: true, sep: ",")
        .map { sample ->
        [sample.name, file(sample.fastq)]}
        .set{ concatenate_ch }

}

// Concatenate reads by sample name
process concatenate {
  container "ubuntu:16.04"
  cpus 1
  memory "4 GB"
  
  input:
  set sample_name, file(fastq_list) from concatenate_ch.groupTuple()
  
  output:
  set sample_name, file("${sample_name}.fastq.gz") into count_reads, filter_host_ch

  afterScript "rm *"

  """
#!/bin/bash

set -e

cat ${fastq_list} > TEMP && mv TEMP ${sample_name}.fastq.gz

  """

}
// Count the number of input reads
process countReads {
  container "ubuntu:16.04"
  cpus 1
  memory "4 GB"
  
  input:
  set sample_name, file(fastq) from count_reads
  
  output:
  file "${sample_name}.countReads.csv" into total_counts

  afterScript "rm *"

  """
#!/bin/bash

set -e

n=\$(gunzip -c "${fastq}" | awk 'NR % 4 == 1' | wc -l)
echo "${sample_name},total_reads,\$n" > "${sample_name}.countReads.csv"

  """

}

if (params.paired || params.interleaved){
  extra_bwa_flag = " -p "
  samtools_filter_unmapped = "-F 2"
  samtools_filter_mapped = "-f 2"
}
else {
  extra_bwa_flag = " "
  samtools_filter_unmapped = "-f 4"
  samtools_filter_mapped = "-F 4"
}

// Filter out any reads that align to the host
process filterHostReads {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 8
  memory "16 GB"

  input:
  file host_genome_tar
  set sample_name, file(fastq) from filter_host_ch
  val min_qual from params.min_qual
  val extra_bwa_flag
  val samtools_filter_unmapped
  val threads from 8
  
  output:
  set sample_name, file("${sample_name}.filtered.fastq.gz") into align_ribo_ch, align_genome_ch, count_nonhuman

  afterScript "rm *"

  """
#!/bin/bash

set -e

# Untar the host genome
tar xvf ${host_genome_tar}

host_genome_name=\$(echo ${host_genome_tar} | sed 's/.tar//')

[[ -s \$host_genome_name ]]

# Align with BWA and save the unmapped BAM
bwa mem -T ${min_qual} -t ${threads}${extra_bwa_flag}\$host_genome_name ${fastq} > ${sample_name}.bam
echo "Number of alignments: \$(samtools view ${sample_name}.bam | wc -l)"

cat ${sample_name}.bam | samtools view ${samtools_filter_unmapped} > ${sample_name}.unmapped.bam
echo "Number of unmapped alignments: \$(samtools view ${sample_name}.unmapped.bam | wc -l)"

samtools fastq ${sample_name}.unmapped.bam | \
gzip -c \
> ${sample_name}.filtered.fastq.gz

    """

}


// Count the number of non-human reads
process countNonhumanReads {
  container "ubuntu:16.04"
  cpus 1
  memory "4 GB"
  
  input:
  set sample_name, file(fastq) from count_nonhuman
  
  output:
  file "${sample_name}.countNonhumanReads.csv" into nonhuman_counts

  afterScript "rm *"

  """
#!/bin/bash

set -e

n=\$(gunzip -c "${fastq}" | awk 'NR % 4 == 1' | wc -l)
echo "${sample_name},nonhuman_reads,\$n" > "${sample_name}.countNonhumanReads.csv"

  """

}


// Align reads against all ribosomes
process alignRibosomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus { 4 * task.attempt }
  memory { 8.GB * task.attempt }
  
  errorStrategy "retry"

  input:
  file ribosome_tar
  set sample_name, file(input_fastq) from align_ribo_ch
  val min_qual from params.min_qual
  val samtools_filter_mapped
  val threads from 4
    
  output:
  set sample_name, file("${sample_name}.ribosome.bam") into ribo_coverage_ch

  afterScript "rm *"

  """
#!/bin/bash

set -e

# Untar the indexed ribosome database
tar xvf ${ribosome_tar}

# Align with BWA and remove unmapped reads
bwa mem -T ${min_qual} -a -t ${threads * task.attempt}${extra_bwa_flag}${params.database_prefix}.ribosomes.fasta ${input_fastq} | samtools view -b ${samtools_filter_mapped} - -o ${sample_name}.ribosome.bam

    """

}

// Calculate the coverage of each ribosome reference
process riboCoverage {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 1
  memory "4 GB"

  input:
  set sample_name, file(bam) from ribo_coverage_ch
  
  output:
  set sample_name, file("${sample_name}.ribosome.pileup"), file("${sample_name}.ribosome.idxstats") into ribo_hits_ch

  afterScript "rm *"

  """
#!/bin/bash

set -e

samtools sort ${bam} > ${bam}.sorted
samtools index ${bam}.sorted
samtools mpileup ${bam}.sorted > ${sample_name}.ribosome.pileup
samtools idxstats ${bam}.sorted > ${sample_name}.ribosome.idxstats
rm ${bam}.sorted ${bam}

  """

}

// Figure out which genomes have sufficient ribosomal coverage
process pickGenomes {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  cpus 1
  memory "4 GB"

  input:
  set sample_name, file(sample_pileup), file(sample_idxstats) from ribo_hits_ch
  file ribosome_tsv
  val min_cov_pct from params.min_cov_pct
  
  output:
  file "${sample_name}.genomes.txt" into sample_genomes_txt

  afterScript "rm *"

  """
#!/usr/bin/env python3
import os
import json
import pandas as pd

# Read in a file with the length of each reference
idxstats = pd.read_csv("${sample_idxstats}", sep="\\t", header=None)
ref_len = idxstats.set_index(0)[1].apply(int)

# Read in a file with a list of the positions covered in the alignment
if os.stat("${sample_pileup}").st_size > 0:
    pileup = pd.read_csv("${sample_pileup}", sep="\\t", header=None)

    # Calculate the coverage as the number of covered bases divided by the length
    ref_cov = pileup.groupby(0).apply(len) / ref_len
    ref_cov = ref_cov.dropna()
    ref_cov.sort_values(ascending=False, inplace=True)
    
    # Find those ribosomes with the minimum threshold met
    detected_ribosomes = ref_cov.index.values[ref_cov >= (float("${min_cov_pct}") / 100)]
    print("\\nDetected ribosomes:")
    print("\\n".join(detected_ribosomes))

    # Read in a list of which genomes match which ribosomes
    genome_list = pd.read_csv("${ribosome_tsv}", sep="\\t", header=None)
    print(genome_list)

    # Get those genomes containing the detected ribosomes
    detected_genomes = list(set(genome_list.loc[genome_list[1].isin(detected_ribosomes), 0]))
    print("\\nDetected genomes:")
    print("\\n".join(detected_genomes))

else:
    detected_genomes = []

# Write out to a file
with open("${sample_name}.genomes.txt", "wt") as fo:
    fo.write("\\n".join(detected_genomes))
"""

}

// Join together the complete list of genomes that were detected across all samples
process aggregateGenomes {
  container "ubuntu:16.04"
  cpus 1
  memory "4 GB"

  input:
  file txt from sample_genomes_txt.collect()
  
  output:
  file "detected.genomes.txt" into genome_hits_ch, gff_hits_ch

  afterScript "rm *"

"""
#!/bin/bash

set -e

for fp in ${txt}; do cat \$fp; echo; done | sort -u | sed '/^\$/d' > TEMP && mv TEMP detected.genomes.txt
"""

}


// Filter down to the genomes detected for each sample
process filterGenomes {
  container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
  cpus 2
  memory "8 GB"

  input:
  file genome_fasta
  file genome_tsv
  file sample_genomes from genome_hits_ch
  publishDir "${params.output_folder}/ref"
  
  output:
  file "filtered.ref.fasta" into filtered_ref_fasta

  afterScript "rm *"

  """
#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip

# Read in the genomes needed for this sample
sample_genomes = open("${sample_genomes}").readlines()
sample_genomes = [fp.rstrip("\\n") for fp in sample_genomes]

# Figure out which headers that corresponds to
genome_headers = dict()
all_headers = set([])
for line in gzip.open("${genome_tsv}", "rt").readlines():
    line = line.rstrip("\\n").split("\\t")
    if len(line) != 2:
        continue
    genome, header = line
    assert header not in all_headers, "Found duplicate header: " + header
    genome_headers[genome] = genome_headers.get(genome, [])
    genome_headers[genome].append(header)

sample_headers = set([
    header
    for genome in sample_genomes
    for header in genome_headers[genome]
])

# Extract the sequences from the FASTA
n_written = 0
with open("filtered.ref.fasta", "wt") as fo:
    for header, seq in SimpleFastaParser(gzip.open("${genome_fasta}", "rt")):
        header = header.split(" ")[0].split("\\t")[0]
        if header in sample_headers:
            fo.write(">" + header + "\\n" + seq + "\\n")
            n_written += 1
assert n_written == len(sample_headers), (n_written, len(sample_headers))

  """

}


// Filter down to the annotations for the genomes detected for each sample
process filterGFF {
  container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
  cpus 1
  memory "4 GB"
  publishDir "${params.output_folder}/ref"

  input:
  file all_gff
  file genome_tsv
  file sample_genomes from gff_hits_ch
  val min_cov_pct from params.min_cov_pct
  
  output:
  file "filtered.ref.gff.gz" into filtered_gff

  afterScript "rm *"

  """
#!/usr/bin/env python3
import gzip

# Read in the genomes needed for this sample
sample_genomes = open("${sample_genomes}").readlines()
sample_genomes = [fp.rstrip("\\n") for fp in sample_genomes]
assert len(sample_genomes) > 0, "No genomes passed the filtering threshold of %d percent" % ("${min_cov_pct}")

# Figure out which headers that corresponds to
genome_headers = dict()
all_headers = set([])
for line in gzip.open("${genome_tsv}", "rt").readlines():
    line = line.rstrip("\\n").split("\\t")
    if len(line) != 2:
        continue
    genome, header = line
    assert header not in all_headers, "Found duplicate header: " + header
    genome_headers[genome] = genome_headers.get(genome, [])
    genome_headers[genome].append(header)

sample_headers = set([
    header
    for genome in sample_genomes
    for header in genome_headers[genome]
])

# Extract the sequences from the GFF
n_written = 0
with gzip.open("filtered.ref.gff.gz", "wt") as fo:
    with gzip.open("${all_gff}", "rt") as fi:
        for line in fi:
            if line[0] == '#':
                continue
            if line.split("\\t")[0] in sample_headers:
                fo.write(line)
                n_written += 1
assert n_written > 0, "No annotations were found for the filtered genomes"

  """

}

// Make an indexed database of the genomes for each sample
process indexGenomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 4
  memory "8 GB"
  publishDir "${params.output_folder}/ref"

  input:
  file filtered_ref_fasta
  
  output:
  file "filtered.ref.fasta.tar" into ref_fasta_tar
  file "${filtered_ref_fasta}.fai"

  afterScript "rm *"

  """
#!/bin/bash
set -e

# Index the selected genomes
bwa index "${filtered_ref_fasta}"
samtools faidx "${filtered_ref_fasta}"

# Tar up the index
tar cvf filtered.ref.fasta.tar filtered.ref.fasta*
    """

}


// Align reads against selected genomes
process alignGenomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus { 4 * task.attempt }
  memory { 16.GB * task.attempt}
  publishDir "${params.output_folder}/bam"
  errorStrategy "retry"

  input:
  set sample_name, file(input_fastq) from align_genome_ch
  file ref_fasta_tar
  val min_qual from params.min_qual
  val extra_bwa_flag
  val samtools_filter_unmapped
  val threads from 4
  
  output:
  set sample_name, file("${sample_name}.genomes.bam") into count_aligned
  set sample_name, file("${sample_name}.genomes.pileup.gz") into genome_pileup

  afterScript "rm *"

  """
#!/bin/bash
set -e

# Untar the indexed genome database
tar xvf ${ref_fasta_tar}

# Align with BWA and remove unmapped reads
bwa mem -T ${min_qual} -a -t ${threads * task.attempt}${extra_bwa_flag}filtered.ref.fasta ${input_fastq} | samtools view -b ${samtools_filter_mapped} - -o ${sample_name}.genomes.bam

samtools sort ${sample_name}.genomes.bam > ${sample_name}.genomes.bam.sorted
mv ${sample_name}.genomes.bam.sorted ${sample_name}.genomes.bam
samtools index ${sample_name}.genomes.bam
samtools mpileup ${sample_name}.genomes.bam | gzip -c > ${sample_name}.genomes.pileup.gz

echo "Number of aligned bases: \$(gunzip -c ${sample_name}.genomes.pileup.gz | wc -l)"

(( \$(gunzip -c ${sample_name}.genomes.pileup.gz | wc -l) > 0 ))

    """

}

// Count the number of aligned reads
process countAlignedReads {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 4
  memory "8 GB"
  
  input:
  set sample_name, file(bam) from count_aligned
  
  output:
  file "${sample_name}.countMapped.csv" into mapped_counts

  afterScript "rm *"

  """
#!/bin/bash

set -e

n=\$(samtools view "${bam}" | cut -f 1 | sort -u | wc -l)
echo "${sample_name},mapped_reads,\$n" > "${sample_name}.countMapped.csv"

  """

}


// Calculate summary metrics for each sample across all genomes
process summarizeAlignments {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  cpus 36
  memory "60 GB"

  input:
  set sample_name, file(sample_pileup) from genome_pileup
  file filtered_gff
  file genome_table
  
  output:
  file "${sample_name}.summary.csv" into sample_results

  afterScript "rm *"

  """
#!/usr/bin/env python3
import os
import json
import gzip
import pandas as pd

# Read in the pileup
# Calculate the depth per base
print("Reading in the pileup file and calculating base depth")
base_depth = dict([
    (reference, reference_pileup.set_index(1)[3].to_dict())
    for reference, reference_pileup in pd.read_csv(
        "${sample_pileup}", 
        header=None, 
        sep="\\t", 
        usecols=[0, 1, 3]
    ).groupby(0)
])

# Read in the list of organism names for each reference
print("Reading in the organism names")
org_names = pd.read_csv(
    "${genome_table}", 
    sep="\\t", 
    header=None
).set_index(1)[0]

all_references = set(org_names.index.values)

# Read in the GFF annotations
print("Reading in the GFF annotations")
annot = []
for line in gzip.open("${filtered_gff}", "rt"):
    if line[0] == '#':
        continue
    line = line.split("\\t")

    if line[0] not in all_references:
        continue
    
    # Get the gene name
    gene_desc = dict([
        field.split("=", 1)
        for field in line[8].split(";")
        if len(field.split("=", 1)) == 2
    ])

    if "ID" not in gene_desc:
        continue
    
    annot.append(dict([
        ("type", line[2]),
        ("reference", line[0]),
        ("start", int(line[3])),
        ("end", int(line[4])),
        ("ID", gene_desc["ID"]),
        ("name", gene_desc.get("Name", gene_desc.get("name"))),
        ("product", gene_desc.get("product"))
    ]))

# Format as a DataFrame
annot = pd.DataFrame(annot)

# Subset to a few types of features
annot = annot.query(
    "type in ['CDS', 'tRNA', 'ncRNA', 'rRNA']"
)

# Add the organism name
annot["organism"] = annot["reference"].apply(org_names.get)

# Precompute the length of each feature
annot["length"] = 1 + annot["end"] - annot["start"]
assert (annot["length"] > 0).all()

# Compute the depth of sequencing for each feature
print("Computing depth of sequencing per gene")
annot["depth"] = annot.apply(
    lambda r: sum([
        base_depth.get(r["reference"], dict()).get(ix, 0)
        for ix in range(r["start"], r["end"] + 1)
    ])/ (r["length"]),
    axis=1
)

# Add the sample name
annot["sample"] = "${sample_name}"

# Write out to a file
print("Writing out to a file")
annot.to_csv("${sample_name}.summary.csv", sep=",", index=None)
print("Done")
"""

}

// Combine results across all genomes
process finalResults {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  cpus 4
  memory "8 GB"
  publishDir "${params.output_folder}"

  input:
  file "*" from sample_results.collect()
  
  output:
  file "${params.output_prefix}.*.csv"

  afterScript "rm *"

  """
#!/usr/bin/env python3
import os
import pandas as pd

df = pd.concat([
    pd.read_csv(fp)
    for fp in os.listdir(".")
    if fp.endswith(".csv")
])

# Write out the complete set of data
df.to_csv("${params.output_prefix}.all.csv", index=None, sep=",")

# Summarize for each organism
# Calculate the weighted average depth for genes broken out by type
summary_df = []
for ix, sub_df in df.groupby(["organism", "sample"]):
    i = dict([
        ("organism", ix[0]),
        ("sample", ix[1])
    ])
    for t, type_df in sub_df.groupby("type"):
        i[t] = (type_df["depth"] * type_df["length"]).sum() / type_df["length"].sum()
    summary_df.append(i)
summary_df = pd.DataFrame(summary_df).set_index(["organism", "sample"])

# Calculate the ratio of rRNA to CDS for each organism & sample
summary_df["rRNA_CDS_ratio"] = summary_df["rRNA"] / summary_df["CDS"]

for t in summary_df.columns.values:
    summary_df.reset_index().pivot_table(
        index="sample",
        columns="organism",
        values=t
    ).reset_index().to_csv(
        "${params.output_prefix}." + t + ".csv",
        index=None,
        sep=","
    )

# Make a column with the annotation of each gene
df["annotation"] = df.apply(
    lambda r: "%s: %s (%s)" % (r["name"], r["product"], r["ID"]),
    axis=1
)

# For each organism, print out the depth of sequencing across all samples
for org, org_df in df.groupby("organism"):
    org_df.pivot_table(
        index="sample",
        columns="annotation",
        values="depth"
    ).reset_index().to_csv(
        "${params.output_prefix}." + org + ".csv",
        index=None,
        sep=","
    )

"""

}


// Combine results across all genomes
process mappingSummary {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  cpus 4
  memory "8 GB"
  publishDir "${params.output_folder}"

  input:
  file "*" from total_counts.collect()
  file "*" from nonhuman_counts.collect()
  file "*" from mapped_counts.collect()
  
  output:
  file "${params.output_prefix}.mapping_summary.csv"

  afterScript "rm *"

  """
#!/usr/bin/env python3
import os
import pandas as pd

df = pd.concat([
    pd.read_csv(fp, header=None, names=["sample", "metric", "value"])
    for fp in os.listdir(".")
    if fp.endswith(".csv")
])

df = df.pivot_table(
    index="sample",
    columns="metric",
    values="value"
)
df.reset_index().to_csv(
    "${params.output_prefix}.mapping_summary.csv",
    index=None,
    sep=","
)

"""

}

