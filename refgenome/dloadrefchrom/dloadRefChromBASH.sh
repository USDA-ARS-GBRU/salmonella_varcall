#!/bin/bash

# This script uses the program ncbi-acc-download to download a list of Genbank
# FASTA files for S. enterica chromosomes, which will be used to create a
# reference graph.

# Script written by John Christian Gaby on 11/18/2020

date

input="/project/gbru_fy21_salmonella/data/accList/refChromList.txt"

while IFS= read -r line
do
/home/chris.gaby/.local/bin/ncbi-acc-download -v --format fasta "$line"
done < "$input"

date
