#!/bin/bash

# Script written by John Christian Gaby
# Original script written Oct. 4, 2016
# This version was adapted for the Salmonella project on Nov. 24, 2020

# This Unix BASH script makes a distance matrix from MASH distances. Sketches
# must already have been calculated and placed into a directory. This script
# takes MASH sketches of the full-length, contiguous Salmonella enterica 
# chromosomes, which are intended for use in generating a reference genome
# graph, and creates a distance matrix of all pairwise distances between the 
# nearly 600 S. enterica chromosomes, yielding 600 X 600 = 360,000 distances.
# The resulting distance matrix is in a format that may be used as the basis
# for constructing an ordination plot to determine whether there are any
# outlier chromosomes that should not be used for building the reference graph. 


# Smaller, test dataset.
# Specify the same fileset for each variable.
#FILES=/project/gbru_fy21_salmonella/data/mash/sketches/CP0637*.msh
#FILES2=/project/gbru_fy21_salmonella/data/mash/sketches/CP0637*.msh

# Full dataset of nearly 600 sketches.
# Specify the same fileset for each variable.
FILES=/project/gbru_fy21_salmonella/data/mash/sketches/*
FILES2=/project/gbru_fy21_salmonella/data/mash/sketches/*

# Begin with a tab so that the column labels are correctly positioned.
printf '\t' >> salmRefChromMASHDist.tab

# This initial for loop creates the column labels.
for n in $FILES2

  do

  printf '%s' $(/home/chris.gaby/mash-Linux64-v2.2/mash dist $n $n | cut -f 2) >> salmRefChromMASHDist.tab
  printf '\t' >> salmRefChromMASHDist.tab

done

printf '\n' >> salmRefChromMASHDist.tab

# Below is a nested for loop. First, the accession is printed at the first of
# the line followed by a tab. Then, the inner for loop prints the distances 
# for all comparisons of the accession to all other accessions, and thereafter
# the outer for loop proceeds to the next accession. Distances are separated by
# tabs and each line ends with a newline character. 
for f in $FILES

  do

  printf '%s' $(/home/chris.gaby/mash-Linux64-v2.2/mash dist $f $f | cut -f 1) >> salmRefChromMASHDist.tab
  printf '\t' >> salmRefChromMASHDist.tab

  # The inner for loop begins here.
  for n in $FILES2

    do

    printf '%s' $(/home/chris.gaby/mash-Linux64-v2.2/mash dist $f $n | cut -f 3) >> salmRefChromMASHDist.tab
    printf '\t' >> salmRefChromMASHDist.tab

  done

  printf '\n' >> salmRefChromMASHDist.tab

done
