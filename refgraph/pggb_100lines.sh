#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH -A gbru_fy21_salmonella
#SBATCH --time=72:00:00  # walltime limit (HH:MM:SS)
#SBATCH -n 48
#SBATCH --partition=atlas    # standard node(s)
#SBATCH --job-name="salm_pggb"
#SBATCH --output="salm_pggb_%j.out"
#SBATCH --error="salm_pggb_%j.err"


ulimit -s unlimited
export PATH=/home/brian.nadon/local/:/home/brian.nadon/miniconda3/bin:/home/brian.nadon/miniconda3/condabin:/sbin:/usr/sbin:/bin:/usr/bin:/usr/lib64/qt-3.3/bin:/apps/sbin:/apps/bin:/opt/slurm/bin:$PATH
echoerr() { echo "$@" 1>&2; }

#define my files and folders
workingdir=/project/gbru_fy21_salmonella/pg_graph/pggb/100lines/

pggb_img=/project/gbru_fy21_salmonella/software/pggb_latest.simg
nlines=100
lineName=pggb_100lines

fastafile=/project/gbru_fy21_salmonella/pg_graph/pggb/100lines/fasta/concat/100lines.fasta

#do pggb
module load singularity/3.5.2
echoerr "timing pggb creation with pggb -i $fastafile -s 100000 -p 80 -n 20 -t 48"
time singularity exec $pggb_img pggb -i $fastafile -v -S -m -s 100000 -p 80 -n 20 -t 48 -o $workingdir


