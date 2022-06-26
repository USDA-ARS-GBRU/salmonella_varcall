#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -A GBRU
#SBATCH -p atlas
#SBATCH --mail-user=adam.r.rivers@gmail.com
#SBATCH --mail-type=END,FAIL

time ../gatk-4.2.6.1/gatk --java-options "-Xmx300g" GenotypeGVCFs \
  -R GCF_000006945.2.fasta \
  -V gendb://100_gvcf_db \
  -O salm100.vcf.gz \
  --sample-ploidy 1
