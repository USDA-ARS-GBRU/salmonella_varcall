#!/bin/bash

#SBATCH --job-name="dlRefChr"
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --mail-user=jchristiangaby@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"
#SBATCH -A gbru_fy21_salmonella

# This script uses the program ncbi-acc-download to download a list of Genbank
# FASTA files for S. enterica chromosomes, which will be used to create a
# reference graph.

# Script written by John Christian Gaby on 11/18/2020

date

input="/project/gbru_fy21_salmonella/data/accList/refChromListTest10.txt"

while IFS= read -r line
do
/home/chris.gaby/.local/bin/ncbi-acc-download --format fasta "$line"
done < "$input"

date
