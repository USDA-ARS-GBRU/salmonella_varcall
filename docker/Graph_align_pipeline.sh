#!/bin/bash

#Use the following line to set your PATH to include the VG and GraphAligner executables.
export PATH=/path/to/vg:/path/to/GraphAligner:$PATH

#Defining the location of the perl script to compress the pack table:
pack_to_seg=/path/to/pack_table_to_seg_cov.pl

#Parameters time!
#Edit these as needed.
num_threads=2 #Number of CPU threads to use.
filter_threshold=0.98 #the % identity below which to filter out alignments.

#This is a function to allow me to write to stderr to profile process times
#in conjunction with the "time" command
echoerr() { echo "$@" 1>&2; }

#The location of the PGGB 250-line graph. Should be absolute and not relative!
#That is, don't use "." or "..":
#First, define the path where the .gfa and .xg live:
graph_path=/path/to/pggb/graphs
#then, we define the gfa (actual graph):
gfa=${graph_path}/pggb_250.gfa
#now, the .xg index of that graph
#Obtained with vg view -Fv pggb_250.gfa > pggb_250.vg;
#then vg index -x pggb_250.xg pggb_250.vg 
xg=${graph_path}/pggb_250.xg

#The location of the SRA Fastq's. Please note, has to be FASTQ and not .sra!
#If we need to both download and unpack the sra's, let me know, we have to change
#this.
#Replace the X's with the actual SRA accession number.
sra=SRAXXXXXXX
sradir=/path/to/fastqs/
#Set the location of the SRA fastqs for EACH ACCESSION here. I don't know google's 
#SRA configuration so you will certainly have to change this.
fq1=${sradir}${sra}_1.fastq
fq2=${sradir}${sra}_2.fastq
alnstat=${sra}_alignstats.txt
echoerr "FQ: $fq1 , $fq2"
echoerr "graph: $gfa"


#Graph Aligner
#we're enforcing global alignment as discussed. 
echoerr "time GraphAligner -t 40 -g $gfa -f $fq1 -f $fq2 -a ${sra}.gam -x vg > ${sra}_alignstats.txt"
time GraphAligner --global-alignment -t $num_threads -g $gfa -f $fq1 -f $fq2 -a ${sra}.gam -x vg > $alnstat

#Next, filter the resulting GAM alignment at the $filter_threshold identity percent
echoerr "timing filtering"
time vg filter -f -u -r $filter_threshold ${sra}.gam > ${sra}.filtered.gam

#Print out some stats of the unfiltered + filtered alignments 
printf "\n${sra} stats:\n" >> $alnstat
vg stats -a ${sra}.gam >> $alnstat
printf "\n${sra} filtered stats:\n" >> $alnstat
vg stats -a ${sra}.filtered.gam >> $alnstat

#Create the final node coverage table, then compress it using Brian Abernathy's script.
#the vg pack command can be very slow! be patient
echoerr "timing pack table creation, 2 threads, from filtered gam"
time vg pack -t $num_threads -x ${lineName}.xg -g ${sra}.filtered.gam -d > ${sra}.coverage.txt
time perl $pack_to_seg -p ${sra}.coverage.txt | gzip > ${sra}.segcov.gz
