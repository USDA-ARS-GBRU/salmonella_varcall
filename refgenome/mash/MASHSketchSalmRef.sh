#!/bin/bash

# Script written by John Christian Gaby
# 11/24/2020

# This script takes a directory of fasta files of Salmonella enterica
# chromosomes and uses MASH https://mash.readthedocs.io/en/latest/
# to create a reduced, kmer-based representation of the genome 
# called a sketch that may then be used to compute distances between
# genomes or even metagenomes. I our case, we wish to determine whether
# there may be outlier genomes that are either mislabeled, contaminated,
# or low quality to remove them before using the nearly 600 S. enterica
# chromosomes as input for constructing the reference graph.
  

FILES=/project/gbru_fy21_salmonella/data/refChrom/*.fa

for n in $FILES

do

/home/chris.gaby/mash-Linux64-v2.2/mash sketch -s 400 -k 16 -o /project/gbru_fy21_salmonella/data/mash/sketches/"$(basename ${n%%.fa})_sketch" $n

done
