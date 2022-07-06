// salmonella-varcall.nf

// Declare syntax version
nextflow.enable.dsl=2
// Script parameters
params.sralist = "$baseDir/docker/v1/3SentericaSRAAccessions.txt"
params.refgenome = "$baseDir/refvcf/GCF_000006945.2.fasta.gz"
params.refgraph = "$baseDir/refgraph/pggb_100.gfa"
params.refvcf = "$baseDir/refvcf/salm100.vcf.gz"
params.refvcftbi = "$baseDir/refvcf/salm100.vcf.gz.tbi"
params.threads = 8


log.info """\
Salmonella Varcall  -  N F
================================
refgenome   : $params.refgenome
sralist   : $params.sralist
"""

// setup

process bwa_index {
  container 'biocontainers/bwa:v0.7.17_cv1'
  input:
    path refgenome
  output:
    path refgenome

    """
    # index
    bwa index $refgenome
    """
}

process gatk_createdict {
  container "broadinstitute/gatk:4.2.6.1"
  input:

    path refgenome
  output:
    path  "refgenome.dict", emit: refdict

  """
  # write sequence dictioany
  java -jar picard.jar CreateSequenceDictionary -REFERENCE $refgenome  -OUTPUT "refgenome.dict"
  """
}

// some map process is needed to sent each sraid from the params.sraidlist to a worker node.

process fetch_SRA {
  machineType 'n2-standard-8'
  container "ncbi/sra-tools"
  input:
    val sraid
  output:
    path "${sraid}_1.fastq", emit:  forwardreads
    path "${sraid}_2.fastq", emit: reversereads

    """
    fasterq-dump $sraid --threads $params.threads
    """
}

process  map_reads {
  container 'biocontainers/bwa:v0.7.17_cv1'
  input:
    val sraid
    path forwardreads
    path reversereads
    path refgenome
    val threads

  output:
      path "${sraid}.sam", emit: sam
  """
  bwa mem -M -t $threads $refgenome $forwardreads $reversereads > ${sraid}.sam
  """
}

process  sort_and_index {
  container "broadinstitute/gatk:4.2.6.1"
  input:
    path sam
    val sraid

  output:
      path "${sraid}.bam", emit: bam
      path "${sraid}.bam.bai" emit: bamindex
  """
  java -jar picard.jar SortSam \
  -INPUT sam \
  -OUTPUT \
  ${sraid}.bam \
  -VALIDATION_STRINGENCY LENIENT \
  -CREATE_INDEX TRUE \
  -SO coordinate

  """
}

process gatk_haplotypecaller {
  container "broadinstitute/gatk:4.2.6.1"
  input:
    path refvcf
    path refgenome
    path bam
    val sraid
    val threads

  output:
    path  "${sraid}.vcf", emit:  vcf

  """
  gatk  \
    HaplotypeCaller  \
    -R  $refgenome \
    --sample-ploidy 1 \
    -I $bam \
    -O ${sraid}.vcf \
    -alleles $refvcf \
    --native-pair-hmm-threads $threads
  """
}

process vg_giraffe {
  container "quay.io/vgteam/vg"
  input:
    path refgraph  //refgraph/index.giraffe.gbz
    path dist      //refgraph/index.dist
    path min       //refgraph/index.min
    val threads
    val sraid
  output:
    path  "sraid.gam" emit: gam
  """
  #create GAM with variant calls

  vg giraffe -t $threads -M 20 -g $refgraph -d $dist -m $min -f SRR9613650_1.fastq -f SRR9613650_2.fastq > ${sraid}.gam


  """
}
workflow {
    Channel.fromPath(params.sralist) \
        | splitText() \
        | fetch_SRA \

}
