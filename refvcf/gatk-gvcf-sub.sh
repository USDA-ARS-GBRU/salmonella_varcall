#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -A GBRU
#SBATCH -p atlas
#SBATCH --array=0-100

echo "SLURM_ARRAY_TASK_ID: " ${SLURM_ARRAY_TASK_ID}

# Load the file names in a directory into an array
file_array=(mapped/*.rg.bam)

# select the input file from the array using the Task id
input_file=${file_array[$SLURM_ARRAY_TASK_ID]}

bname=`basename $input_file .rg.bam`
time ../gatk-4.2.6.1/gatk  HaplotypeCaller  -R  GCF_000006945.2.fasta --emit-ref-confidence GVCF --pcr-indel-model NONE --sample-ploidy 1 -I $input_file  -O gvcf/${bname}.g.vcf

