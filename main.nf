#!/usr/bin/env nextflow

// Minimum coverage of a ribosomal subunit needed for full genome alignment
params.min_cov_pct = 100

// --batchfile is a CSV with two columns, sample name and FASTQ
Channel.from(file(params.batchfile))
       .splitCsv(header: false, sep: ",")
       .map { job ->
       [job[0], file(job[1])]}
       .into{ filter_host_ch }

// --host_genome is a FASTA
host_genome = file(params.host_genome)

// --genome_list is a CSV with three columns, organism, FASTA, and GFF3
Channel.from(file(params.genome_list))
       .splitCsv(header: false, sep: ",")
       .map { job ->
       [job[0], file(job[1]), file(job[2])]}
       .into{ get_ribosome_ch }


Channel.from(file(params.genome_list))
       .splitCsv(header: false, sep: ",")
       .map { job ->
       [job[0], file(job[1])]}
       .into{ get_genome_ch, get_headers_ch }
       
Channel.from(file(params.genome_list))
       .splitCsv(header: false, sep: ",")
       .map { job ->
       [job[0], file(job[2])]}
       .into{ all_gff_ch }


// Index the host genome
process indexHost {

  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 1
  memory "8 GB"
  errorStrategy 'retry'

  input:
  file host_genome
  
  output:
  set "${host_genome}", file("${host_genome}.tar") into indexed_host
  
  """
  bwa index ${host_genome}
  tar cvf ${host_genome}.tar ${host_genome}*
  """

}


// Filter out any reads that align to the host
process filterHostReads {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 8
  memory "8 GB"
  errorStrategy 'retry'

  input:
  set host_genome_name, file(host_genome_tar) from indexed_host
  each set sample_name, file(fastq) from filter_host_ch
  
  output:
  set sample_name, file("${sample_name}.filtered.fastq.gz") into align_ribo_ch, align_genome_ch

  """
#!/bin/bash

# Untar the host genome
tar xvf ${host_genome_tar}

# Align with BWA and convert to FASTQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
bwa mem -t 8 ${host_genome_name} ${input_fastq} | samtools view -b -F 4 -o ${input_fastq}.\$genome_name.bam

rm -f ${input_fastq} ${host_genome_name}*
    """

}

// Extract the ribosomal sequences from the reference genomes provided
process extractRibosomes {
  container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
  cpus 1
  memory "4 GB"
  errorStrategy 'retry'

  input:
  set organism_name, file(fasta), file(gff3) from get_ribosome_ch
  
  output:
  set file("${fasta}.ribosome.fasta"), file("${fasta}.ribosome.tsv") into ribosome_ch

  """
#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Get the location of all ribosomes from the GFF3

# Extract the sequences from the FASTA

# Write out a TSV linking each FASTA header to the organism

  """

}


// Make an indexed database of all ribosomes
process indexRibosomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 8
  memory "8 GB"
  errorStrategy 'retry'

  input:
  set file(fasta), file(tsv) from ribosome_ch.collect()
  
  output:
  file "ribosomes.tar" into ribosome_tar
  file "ribosomes.tsv" into ribosome_tsv

  """
#!/bin/bash

# Concatenate all TSVs
cat *tsv > ALL_TSV && rm *tsv && mv ALL_TSV ribosomes.tsv

# Concatenate all FASTAs
cat *fasta > ALL_FASTA && rm *fasta && mv ALL_FASTA ribosomes.fasta

# Index the ribosomal sequences
bwa index ribosomes.fasta

# Tar up the index
tar cvf ribosomes.tar ribosomes.fasta*

rm ribosomes.fasta*
    """

}


// Extract the headers for each reference genome
process genomeHeaders {
  container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
  cpus 1
  memory "4 GB"
  errorStrategy 'retry'

  input:
  set organism_name, file(fasta) from get_headers_ch
  
  output:
  file "${organism_name}.headers.tsv" into genome_headers
  file "${organism_name}.filepath" into genome_paths

  """
#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser

def safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode=mode)
    return open(fp, mode=mode)

# Extract the headers from the FASTA and write out to TSV
with open("${organism_name}.headers.tsv", "wt") as fo:
    for header, seq in SimpleFastaParser(safe_open("${fasta}")):
        header = header.split(" ")[0].split("\\t")[0].rstrip("\\n")
        fo.write("${organism_name}\\t" + header + "\\n")

# Write out the file name
with open("${organism_name}.filepath", "wt") as fo:
    fo.write("${fasta}\\n")

  """

}


// Combine all of the headers and FASTAs
process concatGenomes {
  container "ubuntu:16.04"
  cpus 1
  memory "4 GB"
  errorStrategy 'retry'

  input:
  set organism_name, file(fasta) from get_genome_ch.collect()
  file tsv from genome_headers.collect()
  file headers from genome_headers.collect()
  file filepath from genome_paths.collect()
  
  output:
  set file("genomes.fasta"), file("genomes.tsv") into all_genomes
  file "genomes.tsv" into genome_table

  """
#!/bin/bash

set -e

cat *filepath | while read fp; do

    [[ -s "${fp}" ]]

    gzip -t "${fp}" && gunzip -c "${fp}" || cat "${fp}"

done > genomes.fasta

cat *tsv | sed '/^$/d' > genomes.tsv

  """

}


// Combine all of the GFF3 files
process concatGFF {
  container "ubuntu:16.04"
  cpus 1
  memory "4 GB"
  errorStrategy 'retry'

  input:
  file gff from all_gff_ch.collect()
  
  output:
  file "genomes.gff3" into all_gff

  """
#!/bin/bash

set -e

cat *gff3 > genomes.gff3

  """

}


// Align reads against all ribosomes
process alignRibosomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 8
  memory "8 GB"
  errorStrategy 'retry'

  input:
  file ribosome_tar
  each set sample_name, file(input_fastq)
  
  output:
  set sample_name, file("${sample_name}.ribosome.bam") into ribo_coverage_ch

  """
#!/bin/bash

# Untar the indexed ribosome database
tar xvf ${ribosome_tar}

# Align with BWA and remove unmapped reads
bwa mem -a -t 8 ribosome.fasta ${input_fastq} | samtools view -b -F 4 -o ${sample_name}.ribosome.bam

rm -f ${input_fastq} ribosome.fasta*
    """

}

// Calculate the coverage of each ribosome reference
process riboCoverage {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 1
  memory "4 GB"
  errorStrategy 'retry'

  input:
  set sample_name, file(ribosome_bam) from ribo_coverage_ch
  
  output:
  set sample_name, file("${sample_name}.ribosome.pileup"), file("${sample_name}.ribosome.idxstats") into ribo_hits_ch

  """
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
  errorStrategy 'retry'

  input:
  each set sample_name, file(sample_pileup), file(sample_idxstats) from ribo_hits_ch
  file ribosome_tsv
  val min_cov_pct from params.min_cov_pct
  
  output:
  set sample_name, file("${sample_name}.genomes.txt") into genome_hits_ch

  """
#!/usr/bin/env python3
import os
import json
import pandas as pd

# Read in a file with the length of each reference
idxstats = pd.read_csv("${sample_idxstats}", sep="\\t", header=None)
ref_len = idxstats.set_index(0)[1]

# Read in a file with a list of the positions covered in the alignment
pileup = pd.read_csv("${sample_pileup}", sep="\\t", header=None)

# Calculate the coverage as the number of covered bases divided by the length
ref_cov = pileup.groupby(0).shape[0] / ref_len
print(ref_cov.sort_values().tail())

# Find those ribosomes with the minimum threshold met
detected_ribosomes = ref_cov.index.values[ref_cov >= (float("${min_cov_pct}") / 100)]
print("\\nDetected ribosomes:")
print("\\n".join(detected_ribosomes))

# Read in a list of which genomes match which ribosomes
genome_list = pd.read_csv("${ribosome_tsv}", sep="\\t")

# Get those genomes containing the detected ribosomes
detected_genomes = list(set(genome_list.loc[genome_list[1].isin(detected_ribosomes)]))
print("\\nDetected genomes:")
print("\\n".join(detected_genomes))

# Write out to a file
with open("${sample_name}.genomes.txt", "wt") as fo:
    fo.write("\\n".join(detected_genomes))
"""

}


// Filter down to the genomes detected for each sample
process filterGenomes {
  container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
  cpus 1
  memory "4 GB"
  errorStrategy 'retry'

  input:
  set file(genome_fasta), file(genome_tsv) from all_genomes
  each set sample_name, file(sample_genomes) from genome_hits_ch
  
  output:
  set sample_name, file("${sample_name}.ref.fasta") into index_sample_ref_ch

  """
#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Read in the genomes needed for this sample
sample_genomes = open("${sample_genomes}").readlines()
sample_genomes = [fp.rstrip("\\n") for fp in sample_genomes]

# Figure out which headers that corresponds to
genome_headers = dict()
all_headers = set([])
for line in open("${genome_tsv}", "rt").readlines():
    line = line.rstrip("\\n").split("\\t")
    if len(line) != 2:
        continue
    genome, header = line
    assert header not in all_headers, "Found duplicate header: " + header
    genome_headers[genome] = genome_headers.get(genome, [])
    genome_headers[genome].append(header)

sample_headers = set([
    header
    for genme in sample_genomes
    for header in genome_headers[genome]
])

# Extract the sequences from the FASTA
with open("${sample_name}.ref.fasta", "wt") as fo:
    for header, seq in SimpleFastaParser(open("${genome_fasta}")):
        if header in sample_headers:
            fo.write(">" + header + "\\n" + seq + "\\n")

  """

}


// Make an indexed database of the genomes for each sample
process indexGenomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 8
  memory "8 GB"
  errorStrategy 'retry'

  input:
  set sample_name, file(sample_fasta) from index_sample_ref_ch
  
  output:
  set sample_name, file("${sample_name}.ref.fasta.tar") into align_genome_ref_ch

  """
#!/bin/bash

# Index the selected genomes
bwa index "${sample_fasta}"

# Tar up the index
tar cvf ${sample_name}.ref.fasta.tar ${sample_name}.ref.fasta*
    """

}


// Align reads against selected genomes
process alignGenomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 8
  memory "8 GB"
  errorStrategy 'retry'

  input:
  set sample_name, file(input_fastq), file(ref_fasta_tar) from align_genome_ch.join(align_genome_ref_ch)
  
  output:
  file "${sample_name}.genomes.bam" into genome_bam
  set sample_name, file("${sample_name}.genomes.pileup") into genome_pileup

  """
#!/bin/bash

# Untar the indexed genome database
tar xvf ${sample_name}.ref.fasta.tar

# Align with BWA and remove unmapped reads
bwa mem -a -t 8 ${sample_name}.ref.fasta ${input_fastq} | samtools view -b -F 4 -o ${sample_name}.genomes.bam

samtools sort ${sample_name}.genomes.bam > ${sample_name}.genomes.bam.sorted
mv ${sample_name}.genomes.bam.sorted ${sample_name}.genomes.bam
samtools index ${sample_name}.genomes.bam
samtools mpileup ${sample_name}.genomes.bam > ${sample_name}.genomes.pileup

    """

}


// Calculate summary metrics for each sample across all genomes
process summarizeResults {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  cpus 1
  memory "4 GB"
  errorStrategy 'retry'

  input:
  each set sample_name, file(sample_pileup) from genome_pileup
  file genome_table
  
  output:
  file "${sample_name}.summary.json" into sample_results

  """
#!/usr/bin/env python3
import os
import json
import pandas as pd

# SUMMARIZE RESULTS BY GENE, AND BY TYPE OF GENE

"""

}

// Combine results across all genomes
process summarizeResults {
  container "quay.io/fhcrc-microbiome/python-pandas:v0.24.2"
  cpus 1
  memory "4 GB"
  errorStrategy 'retry'
  publishDir "${params.output_folder}"

  input:
  file all_summary_json from sample_results.collect()
  
  output:
  file "${params.output_name}" into sample_results

  """
#!/usr/bin/env python3
import os
import json
import pandas as pd

# COMBINE RESULTS ACROSS ALL SAMPLES

"""

}