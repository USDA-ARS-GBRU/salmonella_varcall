# salmonella_varcall
Development of a variant calling pipeline for *Salmonella* genomic data.

## Folder Contents

This repository contains the following folders:

1. `.github/workflows` This folder contains the `docker-publish.yml` file that specifies the github action for creating the docker image from the dockerfile and other files contained in the repository.
2. `amr` This folder contains the Python script that downloads AntiMicrobial Resistance (AMR) data from *Salmonella* records in NCBI. The folder also contains a comma separated values (CSV) file of the downloaded AMR data.
3. `docker` This folder contains the dockerfile that specifies instructions to create the docker image, which includes all software, files, and dependencies required to execute conventional SNP calling as well as the novel reference graph-based variant calling. Files unique to this project and thus unavailable elsewhere are also included in this folder and are copied into the docker image via commands in the dockerfile. Note that a Github action has been created for this repository that builds the docker image and deposits it in the Github Container Registry (GHCR) at the location https://github.com/usda-ars-gbru/salmonella_varcall/pkgs/container/salmonella_varcall
4. `metadata` Metadata for NCBI *Salmonella* datasets relevant to this project as well as the code used to process and explore the metadata.
   1. `assembly` Metadata and analysis from the NCBI Assembly database.
   2. `sra` Metadata and analysis from the NCBI Sequence Read Archive (SRA).
5. `nextflow` The `main.nf` Nextflow script in this folder is the pipeline of commands that are executed in the cloud in order to process the *Salmonella* genomic datasets from the SRA and output variant calls.
6. `refgenome` The scripts and files used to download and analyze the *S. enterica* reference chromosomes, which were subsequently used as the basis for constructing the reference graph.
   1. `dloadrefchrom` The scripts and list of Genbank chromosome accessions needed to download the reference chromosomes.
   2. `fixorigin` The script for running the program `Circlator` to set all the fasta files for all of the reference chromosomes to begin at their origin of replication.
   3. `kmedoids` The RMarkdown (with file extension `.Rmd`) file and its output for performing kmedoids clustering of the MASH distances of the reference chromosomes. Lists of Genbank accesions for the medoid representatives are given for a range of cluster sizes. These lists were subsequently used to determine how many reference chromosomes to use for constructing the reference graph.
   4. `mash` The scripts for generating a MASH sketch for each reference chromosome and a pairwise genomic distance matrix as well as the output distance matrix itself, which was subsequently used for outlier identification and kmedoids analyses.
   5. `outlierID` Analysis of the MASH genomic distance matrix to identify outlier chromosomes for exclusion from the reference graph.
7. `refgraph` The custom *Salmonella enterica* reference graph and the code for generating it.
8. `snpcalling` The traditional SNP calling pipelines to which the graph-based variant calling pipeline will be compared.
   1. `freebayes` SNP calling with the FreeBayes program.
   2. `mpileup` SNP calling with the samtools mpileup program.
   3. `refgenome` The reference genome of *Salmonella enterica* type strain LT2 used in the SNP calling.
9. `varcalling` The variant calling pipeline that uses a custom *Salmonella enterica* reference graph.

## Software Used in This Project

The following is a list of the software used in this project:

### Circlator
Circlator was used to set the origin of replication for the reference *S. enterica* chromosomes that were used to create the reference graph.

[Circlator on Github](https://github.com/sanger-pathogens/circlator)

### FreeBayes

### GATK

The Genome Analysis ToolKit is software created by the Broad Institute for analysis of genomic data. It is highly capable and complex. While the software was originally developed for analysis of eukaryotic genomes, it has been increasingly used for the analysis of prokaryotic genomes despite consequential differences in genome biology between the domains. In recognition of the need for a pipeline that takes account of these differences, the GATK team has worked on a microbe-specific implementation of the pipeline called [GATK for microbes](https://gatk.broadinstitute.org/hc/en-us/articles/360060004292-Introducing-GATK-for-Microbes). In addition, unlike FreeBayes and Mpileup (Mpileup is part of the SRA toolkit), GATK allows start to finish variant calling on large numbers of genomes (i.e. thousands of genomes) by using the [GVCF output format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format), which permits [consolidation of the gvcf](https://gatk.broadinstitute.org/hc/en-us/articles/360035889971) output from [variant calling on multiple batches of samples](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode).

[GATK is available from the Broad institute](https://gatk.broadinstitute.org/hc/en-us).

### MASH

MASH was used to compare similarity of the reference *S. enterica* chromosomes.

### Samtools

### PGGB

The [PanGenome Graph Builder(PGGB)](https://github.com/pangenome/pggb) is the software that was used to create the pangenome reference graph from the complete, contiguous, circular *S. enterica* chromosomes. The pangenome reference graphs and the script to run PGGB are found in the `refgraph` folder.

### VG

[VG](https://github.com/vgteam/vg) is software for working with variation graphs. It serves as the principle tool that we use for alignment of genome sequence reads to the pangenome reference graph as well as for subsequent variant calling.

Note that Giraffe, which the VG team describes as a tool for "fast, haplotype-based mapping of reads to a graph", is a command within the VG software. It is used for mapping the genome sequencing reads to the pangenome reference graph.

### BWA

### SRA Toolkit

### FastQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a java-based program with graphical user interface (GUI) and command line interface for evaluating the quality of sequencing data. It can produce both html output with graphs that are navigable in a web browser as well as a simple text file with "pass/fail" indications for each of several quality analyses. The text "pass/fail" table for each *S. enterica* SRA dataset may be parsed to determine whether a dataset should be discarded due to anomalous results from a quality analysis.

### Dashing

[Dashing](https://github.com/dnbaker/dashing) is software similar to MASH that allows for fast genomic distance calculation. The software is used to compute a sketch of all *S. enterica* genomic datasets in the SRA. The sketches may be subsequently used to create a matrix of genome distances, which can then be used for clustering to identify outlier or anomalous datasets as well as for generating a list of cluster representatives as a form of dereplication to reduce the number of highly similar genomes in the collection.
