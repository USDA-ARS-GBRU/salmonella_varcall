#!/bin/bash

# Traditional SNP and variant calling for S. enterica SRA datasets
# against the S. enterica LT2 reference genome using the FreeBayes
# variant detection program

# written by John Christian Gaby
# August 4th, 2021

#wget https://sra-download.ncbi.nlm.nih.gov/traces/sra79/SRZ/014464/SRR14464606/87712.fastq -O SRR14464606_1.fastq

#wget https://sra-download.ncbi.nlm.nih.gov/traces/sra79/SRZ/014464/SRR14464606/87713.fastq -O SRR14464606_2.fastq

cp /home/chrisgaby/Desktop/dockerbuild/AE006468.2.fasta ./

bwa index AE006468.2.fasta

bwa mem -t 6 /home/chrisgaby/Desktop/dockerbuild/AE006468.2.fasta SRR14464606_1.fastq SRR14464606_2.fastq | samtools sort -o SRR14464606MapToLT2ref.bam -

samtools index SRR14464606MapToLT2ref.bam

samtools flagstat --threads 6 SRR14464606MapToLT2ref.bam > SRR14464606MapToLT2refStat.txt

freebayes -f /home/chrisgaby/Desktop/dockerbuild/AE006468.2.fasta -p 1 SRR14464606MapToLT2ref.bam > SRR14464606MapToLT2ref.vcf

# samtools mpileup is deprecated; the commands should be updated to bcftools mpileup
#/home/programs/samtools_1.12/bin/samtools mpileup -f AE006468.2.fasta -g SRR14464606MapToLT2ref.bam | /home/programs/bcftools_1.12/bin/bcftools call -m - > SRR14464606MapToLT2ref.vcf
