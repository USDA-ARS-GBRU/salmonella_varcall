# salmonella_varcall
Development of a variant calling pipeline for *Salmonella* genomic data.

This repository contains the following folders:

1. `docker` The dockerfile that specifies instructions to create the docker image with all software, files, and dependencies required to execute the traditional SNP and reference graph-based variant calling pipeline.
2. `metadata` Metadata for NCBI *Salmonella* datasets relevant to this project as well as the code used to process and explore the metadata.
   1. `assembly` Metadata and analysis from the NCBI Assembly database.
   2. `sra` Metadata and analysis from the NCBI Sequence Read Archive (SRA).
3. `refgraph` The custom *Salmonella enterica* reference graph and the code for generating it.
4. `snpcalling` The traditional SNP calling pipelines to which the graph-based variant calling pipeline will be compared.
   1. `freebayes` SNP calling with the FreeBayes program.
   2. `mpileup` SNP calling with the samtools mpileup program.
   3. `refgenome` The reference genome of *Salmonella enterica* type strain LT2 used in the SNP calling.
5. `varcalling` The variant calling pipeline that uses a custom *Salmonella enterica* reference graph.
