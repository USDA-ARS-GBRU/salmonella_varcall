#!/bin/#!/usr/bin/env bash

# Scripts to  downlaod and process the 100 reference genome data for GATK gene calling

# Download the genomes
for line in `cat 100Chrom_nuccore.txt`; do efetch -db nuccore -id $line -format fasta > $line.fasta; done

# create reads from the genomes used for VG calling

for file in refgenomes/*
  do
    bname=`basename $file .fasta`
    ~/bbmap/randomreads.sh \
    ref=$file \
    out=simreads/"$bname"_1.fq.gz \
    out2=simreads/"$bname"_2.fq.gz \
    paired=true \
    coverage=50 \
    illuminanames=t
  done

# index the data with bwa version 0.7.12
bwa index GCF_000006945.2.fasta


#2. Align the paired reads to reference genome using bwa mem.
#   Note: Specify the number of threads or processes to use using the -t
#   parameter. The possible number of threads depends on the machine where the command will run.

for file in simreads/1/*
  do
    bname=`basename $file _1.fq.gz`
    time bwa mem -M -t 16 GCF_000006945.2.fasta simreads/1/"$bname"_1.fq.gz simreads/2/"$bname"_2.fq.gz > mapped/"$bname".sam
  done

# sort Bam and create bam and bai files
for file in mapped/*
  do
    bname=`basename $file .sam`
    echo "$bname"
    java -Xmx20G -jar ~/picard.jar SortSam -INPUT mapped/"$bname".sam -OUTPUT mapped/"$bname".sorted.bam -VALIDATION_STRINGENCY LENIENT -CREATE_INDEX TRUE -SO coordinate
  done

#6. Add or replace read groups
for file in mapped/*.sorted.bam
  do
    bname=`basename $file .sorted.bam`
    echo "$bname"
    java -Xmx20G -jar ~/picard.jar AddOrReplaceReadGroups -INPUT mapped/"$bname".sorted.bam -OUTPUT mapped/"$bname".rg.bam -LB lib1 -PL Illumina  -PU 1 --SM $bname -CN GBRU  -VALIDATION_STRINGENCY LENIENT -CREATE_INDEX TRUE -SO coordinate
  done

# index the reference

/home/adam.rivers/samtools-1.15.1/samtools faidx GCF_000006945.2.fasta

# write sequencd dictioany
java -Xmx20G -jar ~/picard.jar CreateSequenceDictionary -REFERENCE GCF_000006945.2.fasta  -OUTPUT GCF_000006945.2.dict

#8a. Variant detection Single sample example
# time ../gatk-4.2.6.1/gatk  HaplotypeCaller  -R  GCF_000006945.2.fasta --emit-ref-confidence GVCF --pcr-indel-model NONE --sample-ploidy 1 -I mapped/CP006053.1.rg.bam  -O CP006053.1.g.vcf
# run this in parallel with SLURM
sbatch gatk-gvcf-sub.sh

#Create mapping file for gvcf
for file in gvcf/*.vcf
  do
    bname=`basename $file .g.vcf`
    echo ${bname}$'\t'${file} >> mapping_of_gcvf.txt
  done

# after gvcf calling combine reads
time ../gatk-4.2.6.1/gatk  GenomicsDBImport \
--genomicsdb-workspace-path 100_gvcf_db \
--sample-name-map mapping_of_gcvf.txt \
-L NC_003197.2 \
-L NC_003277.2


#When all GVCF files have been creates call genoypes using:

#time ../gatk-4.2.6.1/gatk --java-options "-Xmx300g" GenotypeGVCFs \
#  -R GCF_000006945.2.fasta \
#  -V gendb://100_gvcf_db \
#  -O salm100.vcf.gz \
#  --sample-ploidy 1

# run this in parallel with SLURM
sbatch jointgenosub.sh


#8b. Or, for previously ascertained snps, (this will be done on new samples in the docker runs)

#time ../gatk-4.2.6.1/gatk  \
#  HaplotypeCaller  \
#  -R  GCF_000006945.2.fasta \
#  --sample-ploidy 1 \
#  -I mapped/CP006053.1.rg.bam \
#  -O CP006053.1.test.vcf \
#  -alleles salm100.vcf.gz

