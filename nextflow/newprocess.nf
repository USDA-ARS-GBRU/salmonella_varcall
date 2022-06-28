// salmonella-varcall.nf

// Declare syntax version
nextflow.enable.dsl=2
// Script parameters
params.sralist = "$baseDir/docker/v1/3SentericaSRAAccession.txt"
params.refgenome = "$baseDir/refvcf/GCF_000006945.2.fasta.gz"
params.refgraph = "$baseDir/refgraph/pggb_100.gfa"
params.refvcf = "$baseDir/refvcf/salm100.vcf.gz"
params.refvcftbi = "$baseDir/refvcf/salm100.vcf.gz.tbi"
params.threads = 8


log.info """\
Salmonella Varcall  -  N F
================================
refgenome   : $params.genome
reads    : $params.reads
variants : $params.variants
denylist : $params.denylist
results  : $params.results
"""

// setup

process bwa_index {
  container: 'biocontainers/bwa:v0.7.17_cv1'
  input:
    path refgraph
    path refgenome
  output:
    none

    """
    # index
    bwa index $refgraph
    samtools faidx $refgenome

    """
}




// some map process is needed to sent each sraid from the params.sraidlist to a worker node.

process fetch_SRA {
  container: "ncbi/sra-tools"
  input:
    path sraid
    path db
  output:
    path "${sraid}_1.fastq"
    path "${sraid}_2.fastq"

    """
    fasterq-dump --threads $params.threads $sraid
    """
}

process  map_reads {
  container: 'biocontainers/bwa:v0.7.17_cv1'
  input:
    path read1
    path read2
    path ref
    path $mapped

  output:
      path "${read1}.sam"
  """
  bwa mem -M -t $params.threads $ref $read1 $read2 > ${read1}.sam
  """
}

process  sort_and_index {
  container: "broadinstitute/gatk:4.2.6.1"
  input:
    path sam

  output:
      path "${read1}.bam"
      path "${read1}.bam.bai"
  """
  java -jar picard.jar SortSam \
  -INPUT sam \
  -OUTPUT \
  mapped/"$bname".sorted.bam \
  -VALIDATION_STRINGENCY LENIENT \
  -CREATE_INDEX TRUE \
  -SO coordinate


  """
}

process gatk_haplotypecaller {
  container: "broadinstitute/gatk:4.2.6.1"
  input:
    path refvcf
    path refgenome
    path bam
    var

  output:
    path  "called.vcf"

  """
  # write sequence dictioany
  java -jar picard.jar CreateSequenceDictionary -REFERENCE $refgraph  -OUTPUT "refgenome.dict"

  gatk  \
    HaplotypeCaller  \
    -R  $refgenome \
    --sample-ploidy 1 \
    -I $bam \
    -O called.vcf \
    -alleles $refvcf
  """
}

process vg_giraffe {
  container: "quay.io/vgteam/vg"
  input:
    path refgraph
    path refgenome
    path bam
    var
  output:
  """
  #create GAM with variant calls
  vg giraffe -t 6 -M 20 -g ${idxbase}.gg -H ${idxbase}.gbwt -m ${idxbase}.min -d ${idxbase}.dist -f ${sradir}${fq1} -i > $gam
  #create the GAF for processing
  vg convert -G $gam $xg > $gaf

  """
}
workflow {
    Channel.fromPath(params.sralist) \
        | splitText() \
        | fetch_SRA \

}
