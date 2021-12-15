#!/bin/bash

# SNP calling for S. enterica SRA datasets against the S. enterica LT2 reference genome
# written by John Christian Gaby
# June 28th, 2021

# run this script in the Docker container (chrisgaby/vcall) working directory, which is /home/SentericaVariantPipeline/

/root/miniconda3/bin/bwa index AE006468.2.fasta
/root/miniconda3/bin/bwa mem AE006468.2.fasta SRR14464606_1.fastq SRR14464606_2.fastq | /home/programs/samtools_1.12/bin/samtools sort -o AE006468ChromMapLT2ref.bam -

/home/programs/samtools_1.12/bin/samtools index AE006468ChromMapLT2ref.bam
/home/programs/samtools_1.12/bin/samtools flagstat AE006468ChromMapLT2ref.bam > AE006468ChromMapLT2refStat.txt

# samtools mpileup is deprecated; the commands should be updated to bcftools mpileup
/home/programs/samtools_1.12/bin/samtools mpileup -f AE006468.2.fasta -g AE006468ChromMapLT2ref.bam | /home/programs/bcftools_1.12/bin/bcftools call -m - > AE006468ChromMapLT2ref.vcf
