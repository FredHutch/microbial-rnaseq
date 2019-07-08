#!/usr/bin/env nextflow

// --database_prefix is the name for the database files
params.database_prefix = "microbial_genomes"

// --output_folder is the location for the database output
params.output_folder = "./"

// --host_genome is a FASTA
host_genome = file(params.host_genome)

// --genome_list is a CSV with three columns, organism, FASTA, and GFF3
Channel.from(file(params.genome_list))
       .splitCsv(header: false, sep: ",")
       .map { job ->
       [job[0], job[1]]}
       .set{ download_genome_ch }

Channel.from(file(params.genome_list))
       .splitCsv(header: false, sep: ",")
       .map { job ->
       [job[0], job[2]]}
       .set{ download_gff_ch }


process downloadGenome {
  container "quay.io/fhcrc-microbiome/wget@sha256:98b90e8bb8a171182505f1e255b0bd85cbbda68f08c08b4877c3fc48e63ac82f"
  cpus 1
  memory "4 GB"
  errorStrategy 'retry'
  
  input:
  set organism_name, fasta_url from download_genome_ch
  
  output:
  set organism_name, file("${organism_name}.fasta.gz") into get_headers_ch, get_ribosome_fasta_ch
  file "${organism_name}.fasta.gz" into get_genome_ch

  afterScript "rm *"

  """
#!/bin/bash

set -e

wget -O ${organism_name}.fasta.gz ${fasta_url}

gzip -t ${organism_name}.fasta.gz

  """
}

process downloadGFF {
  container "quay.io/fhcrc-microbiome/wget@sha256:98b90e8bb8a171182505f1e255b0bd85cbbda68f08c08b4877c3fc48e63ac82f"
  cpus 1
  memory "4 GB"
  errorStrategy 'retry'
  
  input:
  set organism_name, gff_url from download_gff_ch
  
  output:
  set organism_name, file("${organism_name}.gff.gz") into get_ribosome_gff_ch, all_gff_ch

  afterScript "rm *"

  """
#!/bin/bash

set -e

wget -O ${organism_name}.gff.gz ${gff_url}

gzip -t ${organism_name}.gff.gz

  """
}


// Index the host genome
process indexHost {

  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 8
  memory "60 GB"
  publishDir "${params.output_folder}"

  input:
  file host_genome
  
  output:
  file "${host_genome}.tar"
  
  afterScript "rm *"

  """
#!/bin/bash

set -e

bwa index ${host_genome}
tar cvf ${host_genome}.tar ${host_genome}*
  """

}


// Extract the ribosomal sequences from the reference genomes provided
process extractRibosomes {
  container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
  cpus 1
  memory "4 GB"

  input:
  set organism_name, file(fasta), file(gff3) from get_ribosome_fasta_ch.join(get_ribosome_gff_ch)
  
  output:
  set file("${fasta}.ribosome.fasta"), file("${fasta}.ribosome.tsv") into ribosome_ch

  afterScript "rm *"

  """
#!/usr/bin/env python3
import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq

def safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode=mode)
    return open(fp, mode=mode)

# Get the location of all ribosomes from the GFF3
ribosomes = []
for line in safe_open("${gff3}"):
    if line[0] == '#':
        continue
    line = line.split("\\t")
    if line[2] == "rRNA":
        # Get the gene name
        gene_desc = dict([
            field.split("=", 1)
            for field in line[8].split(";")
        ])
        assert "ID" in gene_desc
        # Header, start, end, strand
        assert line[6] in ["+", "-"], line[6]
        assert int(line[3]) < int(line[4])
        ribosomes.append((line[0], int(line[3]), int(line[4]), line[6], gene_desc["ID"]))

# Make sure that the gene names are all unique
n_unique = len(set([f[4] for f in ribosomes]))
assert n_unique == len(ribosomes), (n_unique, len(ribosomes))

# Extract the sequences from the FASTA
# Write out a TSV linking each FASTA header to the organism
n_written = 0
with open("${fasta}.ribosome.fasta", "wt") as fasta_out, open("${fasta}.ribosome.tsv", "wt") as tsv_out:
    for h, s in SimpleFastaParser(safe_open("${fasta}")):
        h = h.split(" ")[0].split("\\t")[0]
        for header, start, end, strand, gene_id in ribosomes:
            if header == h:
                
                gene_sequence = s[start - 1: end]
                if strand == "-":
                    gene_sequence = str(Seq(gene_sequence).reverse_complement())

                gene_name = header + "_" + gene_id
                fasta_out.write(">" + gene_name + "\\n" + gene_sequence + "\\n")
                tsv_out.write("${organism_name}\\t" + gene_name + "\\n")
                n_written += 1

assert n_written == len(ribosomes), (n_written, len(ribosomes))

  """

}


// Make an indexed database of all ribosomes
process indexRibosomes {
  container "quay.io/fhcrc-microbiome/bwa@sha256:2fc9c6c38521b04020a1e148ba042a2fccf8de6affebc530fbdd45abc14bf9e6"
  cpus 4
  memory "8 GB"
  publishDir "${params.output_folder}"

  input:
  file "*" from ribosome_ch.collect()
  val database_prefix from params.database_prefix
  
  output:
  file "${database_prefix}.ribosomes.tar"
  file "${database_prefix}.ribosomes.tsv"

  afterScript "rm *"

  """
#!/bin/bash

set -e

# Concatenate all TSVs
cat *tsv > ALL_TSV && rm *tsv && mv ALL_TSV ${database_prefix}.ribosomes.tsv

# Concatenate all FASTAs
cat *fasta > ALL_FASTA && rm *fasta && mv ALL_FASTA ${database_prefix}.ribosomes.fasta

# Index the ribosomal sequences
bwa index ${database_prefix}.ribosomes.fasta

# Tar up the index
tar cvf ${database_prefix}.ribosomes.tar ${database_prefix}.ribosomes.fasta*
    """

}


// Extract the headers for each reference genome
process genomeHeaders {
  container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
  cpus 1
  memory "4 GB"

  input:
  set organism_name, file(fasta) from get_headers_ch
  
  output:
  file "${organism_name}.headers.tsv.gz" into genome_headers

  afterScript "rm *"

  """
#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip

def safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode=mode)
    return open(fp, mode=mode)

# Extract the headers from the FASTA and write out to TSV
with gzip.open("${organism_name}.headers.tsv.gz", "wt") as fo:
    for header, seq in SimpleFastaParser(safe_open("${fasta}")):
        header = header.split(" ")[0].split("\\t")[0].rstrip("\\n")
        fo.write("${organism_name}\\t" + header + "\\n")

  """

}



process concatGenomes {
  container "ubuntu:16.04"
  cpus 8
  memory "60 GB"
  publishDir "${params.output_folder}"
  
  input:
  file "*" from get_genome_ch.collect()
  val database_prefix from params.database_prefix
  
  output:
  file "${database_prefix}.fasta.gz"
  
  afterScript "rm *"

  """
#!/bin/bash

set -e

cat *fasta.gz >> ${database_prefix}.fasta.gz

  """

}


// Combine all of the headers
process concatHeaders {
  container "ubuntu:16.04"
  cpus 8
  memory "60 GB"
  publishDir "${params.output_folder}"
  
  input:
  file "*" from genome_headers.collect()
  val database_prefix from params.database_prefix
  
  output:
  file "${database_prefix}.tsv.gz"

  afterScript "rm *"

  """

#!/bin/bash

set -e

cat *tsv.gz > ${database_prefix}.tsv.gz

  """

}


// Combine all of the GFF3 files
process concatGFF {
  container "ubuntu:16.04"
  cpus 8
  memory "60 GB"
  publishDir "${params.output_folder}"
  
  input:
  file "*" from all_gff_ch.collect()
  val database_prefix from params.database_prefix
  
  output:
  file "${database_prefix}.gff.gz"

  afterScript "rm *"

  """
#!/bin/bash

set -e

cat *gff.gz > ${database_prefix}.gff.gz

  """

}
