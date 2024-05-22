#!/usr/bin/env bash

# Create new directory for project and change to the directory
mkdir chipseq && cd "$_"

# Download files
fastq-dump --split-files SRR28097737
fastq-dump --split-files SRR28097738
fastq-dump --split-files SRR28097739
fastq-dump --split-files SRR28097740

# Make directories for Quality Control, Mapping and Peaking Calling
mkdir QC_Reports Mapping Peak_Calling

# Quality control using fastqc
fastqc  SRR28097737_1.fastq SRR28097738_1.fastq  -o QC_Reports
fastqc  SRR28097737_1.fastq SRR28097738_1.fastq  -o QC_Reports
fastqc  SRR28097737_1.fastq SRR28097738_1.fastq  -o QC_Reports
fastqc  SRR28097737_1.fastq SRR28097738_1.fastq  -o QC_Reports

# Summarize results with multiqc
multiqc QC_Reports

# Trimming using sickle
sickle se -f SRR28097737_1.fastq -t sanger -o trimmed_SRR28097737_1.fastq -q 20 -l 40
sickle se -f SRR28097738_1.fastq -t sanger -o trimmed_SRR28097738_1.fastq -q 20 -l 40
sickle se -f SRR28097739_1.fastq -t sanger -o trimmed_SRR28097739_1.fastq -q 20 -l 40
sickle se -f SRR28097740_1.fastq -t sanger -o trimmed_SRR28097740_1.fastq -q 20 -l 40

# Download the reference genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.p14.genome.fa.gz

# Index the reference genome
bwa index GRCh38.p14.genome.fa

# Mapping
bwa mem -t 32  GRCh38.p14.genome.fa SRR28097737_1.fastq > Mapping/SRR28097737.sam
bwa mem -t 32  GRCh38.p14.genome.fa SRR28097738_1.fastq > Mapping/SRR28097738.sam
bwa mem -t 32  GRCh38.p14.genome.fa SRR28097739_1.fastq > Mapping/SRR28097739.sam
bwa mem -t 32  GRCh38.p14.genome.fa SRR28097740_1.fastq > Mapping/SRR28097740.sam

# Normalization with samtools

# Convert to sam files

# Sort the data

# Peak calling with mac2
