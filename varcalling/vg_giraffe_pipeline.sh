###################
### USAGE NOTES ###
###################

#USAGE: 
#vg_giraffe_pipeline.sh [name_of_ref_graph] [fastq_1] [fastq_2]

#Any line with [tk] needs to be filled in with the appropriate hard-coded
#directory on the final server/compute node.


gname=$1
fq1=$2
fq2=$3

echoerr() { echo "$@" 1>&2; }

############################
### SET FILE DEFINITIONS ###
############################

#The directory you're running this script in
workingdir=./

#The directory the graph(s) and indices are stored in 
#NB THEY MUST ALL BE NAMED BY $gname! e.g. "pggb_100"
#A reminder: ths should include a .GBWT, .gg, .min, and .dist
graphdir=[tk]

#The location of the script to redistribute duplicate reads
#(Needed because giraffe doesn't handle this by default for some reason)
#https://github.com/USDA-ARS-GBRU/PanPipes/blob/main/scripts/gafRedistribute.pl
dis_script=[tk]

#location of the segcov.pl script 
#https://github.com/brianabernathy/gfa_var_genotyper/blob/main/pack_table_to_seg_cov.pl
segcov_script=[TK]


samplename=${fq1%_1.fastq}
outdir=${workingdir}/${gname}_samplename/

#####################
### BUILD INDICES ###
#####################

#This section is not necessary since we've already built the indices.
#But kept here for recordkeeping.

#cd $graphdir
#vg mod -X 1000 ${gname}.gfa > ${gname}.mod.gfa
#vg gbwt -G ${gname}.mod.gfa -p -o ${gname}.mod.gbwt -g ${gname}.mod.gg
#vg minimizer -g ${gname}.mod.gbwt -i ${gname}.mod.min -G ${gname}.mod.gg
#vg convert -x -g ${gname}.mod.gfa > ${gname}.mod.xg
#vg snarls -T ${gname}.mod.xg > ${gname}.mod.snarls
#vg index -s ${gname}.mod.snarls -j ${gname}.mod.dist ${gname}.mod.xg


#Move to the output dir
mkdir -p $outdir
cd $outdir


########################
### BEGIN ALIGNMENT ####
########################

base=$samplename

#These are dependent on that graph name. Make sure they all match!
#e.g. all graphs/indices are named "pggb_100.mod.gfa" and "pggb_100.mod.gbwt"
#They need to all be in $graphdir as well.
gfa=${graphdir}/${gname}.mod.gfa
idxbase=${graphdir}/${gname}.mod
xg=${idxbase}.xg

out_base="${gname}_${base}_giraffe"
alnstat=${out_base}_alignstats.txt
echoerr "FQ: $fq1 $fq2"
echoerr "graph: $idxbase"

#Setting output names
gam=${out_base}.gam
gaf=${out_base}.gaf

echoerr "timing giraffe, 48 threads, -M 20"
#NOTE: add "-o gaf" for final pipeline maybe? We have to go gaf -> redistribute -> GAM
time vg giraffe -t 48 -M 20 -g ${idxbase}.gg -H ${idxbase}.gbwt -m ${idxbase}.min -d ${idxbase}.dist -f $fq1 -f $fq2 > $gam

#create the GAF for processing
vg convert -G $gam $xg > $gaf

#print out some stats
printf "${out_base} stats:\n" >> $alnstat
vg stats -a $gam >> $alnstat



#randomly distribute duplicate reads
#(reads that align equally to different parts of the genome)
#Note this permanently deletes the original non-fixed GAF and GAM.
echoerr "timing redis"
time perl $dis_script $gaf > tmp.gaf
mv tmp.gaf $gaf
vg convert -F $gaf $xg > $gam

#######################
### VARIANT CALLING ###
#######################
echoerr "timing vg pack"
time vg pack -t 40 -x $xg -g $gam -o ${out_base}.pack -d > ${out_base}.packtable
echoerr "timing vg call"
time vg call -t 40 -p $refname -d 1 -r ${idxbase}.snarls -k ${out_base}.pack -a -s ${out_base} $xg > ${out_base}.vcf
time vg call -t 40 -d 1 -r ${idxbase}.snarls -k ${out_base}.pack -a -s ${out_base} $xg > ${out_base}_allpaths.vcf

#Convert the pack table to our final desired format
perl $segcov_script -p ${out_base}.packtable > ${out_base}.seg.cov