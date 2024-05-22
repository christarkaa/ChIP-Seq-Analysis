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

# Normalisation
# Filter out mitochondrial and unassigned reads
sed '/chrM/d;/random/d;/chrUn/d'  Mapping/SRR28097737.sam >  Mapping/filtered.SRR28097737.sam
sed '/chrM/d;/random/d;/chrUn/d'  Mapping/SRR28097738.sam >  Mapping/filtered.SRR28097738.sam
sed '/chrM/d;/random/d;/chrUn/d'  Mapping/SRR28097739.sam >  Mapping/filtered.SRR28097739.sam
sed '/chrM/d;/random/d;/chrUn/d'  Mapping/SRR28097740.sam >  Mapping/filtered.SRR28097740.sam

# Convert sam to bam files
samtools view -@ 20 -S -b Mapping/filtered.SRR28097737.sam > Mapping/SRR28097737.bam
samtools view -@ 20 -S -b Mapping/filtered.SRR28097738.sam > Mapping/SRR28097738.bam
samtools view -@ 20 -S -b Mapping/filtered.SRR28097739.sam > Mapping/SRR28097739.bam
samtools view -@ 20 -S -b Mapping/filtered.SRR28097740.sam > Mapping/SRR28097740.bam

# Sort the bam files
samtools sort -@ 32 -o Mapping/SRR28097737.sorted.bam Mapping/SRR28097737.bam
samtools sort -@ 32 -o Mapping/SRR28097738.sorted.bam Mapping/SRR28097738.bam
samtools sort -@ 32 -o Mapping/SRR28097739.sorted.bam Mapping/SRR28097738.bam
samtools sort -@ 32 -o Mapping/SRR28097740.sorted.bam Mapping/SRR28097738.bam

# Count number of alignments
samtools view -c Mapping/SRR28097737.sorted.bam
samtools view -c Mapping/SRR28097738.sorted.bam
samtools view -c Mapping/SRR28097739.sorted.bam
samtools view -c Mapping/SRR28097740.sorted.bam

# We won't be doing further normalisation since the alignment reads in the control and treatment files are almost the same

# Rename the files
mv Mapping/SRR28097737.sorted.bam H128-0-Input.bam
mv Mapping/SRR28097738.sorted.bam H128-0-H3K27ac.bam
mv Mapping/SRR28097739.sorted.bam H128-DN5-Input.bam
mv Mapping/SRR28097740.sorted.bam H128-DN5-H3K27ac.bam

# Index the files before peak calling
samtools index H128-0-H3K27ac.bam
samtools index H128-0-Input.bam
samtools index H128-DN5-H3K27ac.bam
samtools index H128-DN5-Input.bam

# Peak calling with macs2
macs2 callpeak -t H128-0-H3K27ac.bam -c H128-0-Input.bam -n Peak_Calling/H128-0-H3K27ac -g hs --bdg -q 0.05 -f BAM
macs2 callpeak -t H128-DN5-H3K27ac.bam -c H128-DN5-Input.bam -n Peak_Calling/H128-DN5-H3K27ac -g hs --bdg -q 0.05 -f BAM
