# salmonella_varcall
Development of a variant calling pipeline for *Salmonella* genomic data.

This repository contains the following folders:

1. `docker` The dockerfile that specifies instructions to create the docker image with all software, files, and dependencies required to execute the traditional SNP and reference graph-based variant calling pipeline.
2. `metadata` Metadata for NCBI *Salmonella* datasets relevant to this project as well as the code used to process and explore the metadata.
   1. `assembly` Metadata and analysis from the NCBI Assembly database.
   2. `sra` Metadata and analysis from the NCBI Sequence Read Archive (SRA).
3. `refgenome` The scripts and files used to download and analyze the *S. enterica* reference chromosomes, which were subsequently used as the basis for constructing the reference graph.
   1. `dloadrefchrom` The scripts and list of Genbank chromosome accessions needed to download the reference chromosomes. 
   2. `fixorigin` The script for running the program `Circlator` to set all the fasta files for all of the reference chromosomes to begin at their origin of replication.
   3. `kmedoids` The Rmd file and its output for performing kmedoids clustering of the MASH distances of the reference chromosomes. Lists of medoid representatives are given for a range of cluster sizes. These lists were subsequently used to determine how many reference chromosomes to use for constructing the reference graph.
   4. `mash` The scripts for generating a MASH sketch for each reference chromosome and a pairwise genomic distance matrix as well as the output distance matrix itself, which was subsequently used for outlier identification and kmedoids analyses.
   5. `outlierID` Analysis of the MASH genomic distance matrix to identify outlier chromosomes for exclusion from the reference graph. 
4. `refgraph` The custom *Salmonella enterica* reference graph and the code for generating it.
5. `snpcalling` The traditional SNP calling pipelines to which the graph-based variant calling pipeline will be compared.
   1. `freebayes` SNP calling with the FreeBayes program.
   2. `mpileup` SNP calling with the samtools mpileup program.
   3. `refgenome` The reference genome of *Salmonella enterica* type strain LT2 used in the SNP calling.
6. `varcalling` The variant calling pipeline that uses a custom *Salmonella enterica* reference graph.
