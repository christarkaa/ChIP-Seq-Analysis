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

# Download the reference genome index for bowtie2
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.p14.genome.fa.gz

# Build the reference genome index
bowtie2-build GRCh38.p14.genome.fa hg38

# Mapping with bowtie 2
bowtie2 -p 10 -x hg38 -U trimmed_SRR28097737_1.fastq -S Mapping/SRR28097737.sam
bowtie2 -p 10 -x hg38 -U trimmed_SRR28097738_1.fastq -S Mapping/SRR28097738.sam
bowtie2 -p 10 -x hg38 -U trimmed_SRR28097739_1.fastq -S Mapping/SRR28097739.sam
bowtie2 -p 10 -x hg38 -U trimmed_SRR28097740_1.fastq -S Mapping/SRR28097740.sam

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
mv Mapping/SRR28097737.sorted.bam Peak_Calling/H128-0-Input.bam
mv Mapping/SRR28097738.sorted.bam Peak_Calling/H128-0-H3K27ac.bam
mv Mapping/SRR28097739.sorted.bam Peak_Calling/H128-DN5-Input.bam
mv Mapping/SRR28097740.sorted.bam Peak_Calling/H128-DN5-H3K27ac.bam

# Index the files before peak calling
samtools index Peak_Calling/H128-0-H3K27ac.bam
samtools index Peak_Calling/H128-0-Input.bam
samtools index Peak_Calling/H128-DN5-H3K27ac.bam
samtools index Peak_Calling/H128-DN5-Input.bam

# Peak calling with macs2
macs2 callpeak -t Peak_Calling/H128-0-H3K27ac.bam -c Peak_Calling/H128-0-Input.bam -n Peak_Calling/H128-0-H3K27ac -g hs --bdg -q 0.05 -f BAM
macs2 callpeak -t Peak_Calling/H128-DN5-H3K27ac.bam -c Peak_Calling/H128-DN5-Input.bam -n Peak_Calling/H128-DN5-H3K27ac -g hs --bdg -q 0.05 -f BAM

# Get chrom.sizes using samtools
fetchChromSizes hg38 > hg38.chromSizes

# Remove the lines in the bedgraph files that are not chromosome names
# View the first column
awk '{print $1}' Peak_Calling/H128-0-H3K27ac_control_lambda.bdg | head
# Preview lines to be removed
awk '$1 ~ /\.[123]$/' Peak_Calling/H128-0-H3K27ac_control_lambda.bdg

# Remove the lines
awk '$1 !~ /\.[123]$/' Peak_Calling/H128-0-H3K27ac_control_lambda.bdg > Peak_Calling/filtered.H128-0-H3K27ac_control_lambda.bdg
awk '$1 !~ /\.[123]$/' Peak_Calling/H128-0-H3K27ac_treat_pileup.bdg > Peak_Calling/filtered.H128-0-H3K27ac_treat_pileup.bdg
awk '$1 !~ /\.[123]$/' Peak_Calling/H128-DN5-H3K27ac_control_lambda.bdg > Peak_Calling/filtered.H128-DN5-H3K27ac_control_lambda.bdg
awk '$1 !~ /\.[123]$/' Peak_Calling/H128-DN5-H3K27ac_treat_pileup.bdg > Peak_Calling/filtered.H128-DN5-H3K27ac_treat_pileup.bdg

#Verify
head Peak_Calling/filtered.H128-0-H3K27ac_control_lambda.bdg

# Convert bedGraph to bigWig
bedGraphToBigWig Peak_Calling/filtered.H128-0-H3K27ac_control_lambda.bdg hg38.chromSizes Peak_Calling/H128-0-H3K27ac_control.bw
bedGraphToBigWig Peak_Calling/filtered.H128-0-H3K27ac_treat_pileup.bdg hg38.chromSizes Peak_Calling/H128-0-H3K27ac_treat.bw
bedGraphToBigWig Peak_Calling/filtered.H128-DN5-H3K27ac_control_lambda.bdg hg38.chromSizes Peak_Calling/H128-DN5-H3K27ac_control.bw
bedGraphToBigWig Peak_Calling/filtered.H128-DN5-H3K27ac_treat_pileup.bdg hg38.chromSizes Peak_Calling/H128-DN5-H3K27ac_trea.bw
