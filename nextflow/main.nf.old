params.index = "/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/1SentericaSRAAccession.txt"

// /home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/10SentericaSRAAccessions.txt

params.project = "SRR12660755"

params.resultdir = '/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/results'

params.outputdir = "/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/work"

params.qualevaldir = "./FastQC2"

LT2ref = Channel.fromPath( './AE006468.2.fasta' )

// projectSRId = params.project

int threads = Runtime.getRuntime().availableProcessors()

process getSRA {
		
	publishDir './SRAs/', mode: 'copy'

	cpus 1

	input:
	path 'accIDs' from params.index
	// val projectID from projectSRId
	
	output:
	file 'sra.txt' into sraIDs
	
	script:
	"""
	tail accIDs -n +2 > sra.txt
	"""
}

sraIDs.splitText().map { it -> it.trim() }.set { singleSRAId }
// singleSRAId.into { singleSRAIdFQD; singleSRAIdMR }

process fastqDump {

	publishDir './SRAs/', mode: 'copy', overwrite: true

	input:
	val id from singleSRAId

	output:
	file '*.fastq.gz' into reads

	//Note that the number of spot IDs should be determined based on read length for a dataset in order to achieve 60X genome coverage

	script:
	"""
	/home/SentericaVariantPipeline/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --clip --skip-technical --maxSpotId 600000 --defline-seq '@\$sn/\$ri' --defline-qual '+\$sn/\$ri' --split-spot --gzip ${id}
	"""
	//600,000 reads needed at 250 X 2 for 60X coverage
	// note that if only one file is output, then in order to separate deflines with "/1" and "/2" to distinguish forward and reverse reads, one must include the "--split-spot" parameter
	//fastq-dump --clip --skip-technical --defline-seq '@$sn/$ri' --defline-qual '+$sn/$ri' --split-files SRR13268953
	// --readids
	// --read-filter pass
	// --split-files
	//parallel-fastq-dump --sra-id $id --threads ${task.cpus} --gzip
	// --maxSpotId 2000000
}

reads.into { reads_dashing; reads_fastqc; reads_gunzip; }

process dashingSketch {

	publishDir './dashing2/', mode: 'move', overwrite: false

	cpus threads

	input:
	file read from reads_dashing
	
	output:
	file '*.hll'
	
	script:
	readName = read.toString() - ~/(\.fastq\.gz)?$/
	
	"""
	/home/SentericaVariantPipeline/dashing_s256 sketch --nthreads $threads -S 10 -k 16 $read
	"""
}

process fastQC {

	publishDir './FastQC2/', mode: 'move', overwrite: false

	cpus threads

	input:
	file read from reads_fastqc
	
	output:
	file '*'
	
	script:
	
	"""
	/home/SentericaVariantPipeline/FastQC/fastqc -t $threads -f fastq --extract -o ./ $read
	"""
}

process indexRef {

	publishDir './', mode: 'move', overwrite: false

	input:
	file ref from LT2ref

	output:
	file '*'
	//into bwa_index

	script:

	"""
	/root/miniconda3/bin/bwa index $ref
	"""
}

process gunzipReads {

	publishDir './SRAs/', mode: 'copy', overwrite: true

	//echo true

	output:
	file '*.fastq' into reads_unzipped

	input:
	file read from reads_gunzip

	script:

	"""
	gunzip -f ${read}
	"""
	//echo ${read}
	//echo ${read.baseName}
}

reads_unzipped.into { reads_unzippedMR; reads_unzippedPGRGvarcall } 

process mapReads {

	publishDir './mappedToRef2/', mode: 'copy', overwrite: false

	cpus threads

	input:
	file read from reads_unzippedMR
	//file '*' from bwa_index.first()
	//val id from singleSRAIdMR

	output:
	file '*.bam' into bamFile

	script:

	"""
	/root/miniconda3/bin/bwa mem /home/SentericaVariantPipeline/AE006468.2.fasta ${read} -t $threads | /home/programs/samtools_1.12/bin/samtools sort -o ${read.simpleName}ChromMappedToLT2ref.bam
	/home/programs/samtools_1.12/bin/samtools index ${read.simpleName}ChromMappedToLT2ref.bam
	"""
}

bamFile.into { bamFile_stats; bamFile_mpileup; bamFile_freebayes }

process mapStats {

	publishDir './mapStats2/', mode: 'move', overwrite: false

	input:
	file bamFile from bamFile_stats

	output:
	file "*.txt"

	script:

	"""
	/home/programs/samtools_1.12/bin/samtools flagstat ${bamFile} > ${bamFile.simpleName}Stats.txt
	"""
}

process mpileup {

	publishDir './mpileup2/', mode: 'move', overwrite: false

	input:
	file bamFile from bamFile_mpileup

	output:
	file "*.vcf"

	script:

	"""
	/home/programs/bcftools_1.12/bin/bcftools mpileup -Ou -f /home/SentericaVariantPipeline/AE006468.2.fasta ${bamFile} | /home/programs/bcftools_1.12/bin/bcftools call -mv -Ov -o ${bamFile.simpleName}_MP.vcf
	"""
	//samtools mpileup -f /home/SentericaVariantPipeline/AE006468.2.fasta -g ${bamFile} | bcftools call -mv - > ./mpileup/${bamFile.simpleName}_MP2.vcf
}

process freebayes {

	publishDir './freebayes2/', mode: 'move', overwrite: false

	input:
	file bamFile from bamFile_freebayes

	output:
	file "*.vcf"

	script:
	"""
	/root/miniconda3/bin/freebayes -f /home/SentericaVariantPipeline/AE006468.2.fasta -p 1 ${bamFile} > ${bamFile.simpleName}_FB2.vcf
	"""
}

process vgGiraffe {

	publishDir './vg2/', mode: 'copy', overwrite: false

	input:
	file read from reads_unzippedPGRGvarcall

	output:
	file '*' into gamfile
	//path "pggb_100_${read.simpleName}.fastq"

	script:

	"""

	/home/programs/vg/vg giraffe -t 6 -M 20 -g /home/SentericaVariantPipeline/100lines/pggb_100.mod.gg -H /home/SentericaVariantPipeline/100lines/pggb_100.mod.gbwt -m /home/SentericaVariantPipeline/100lines/pggb_100.mod.min -d /home/SentericaVariantPipeline/100lines/pggb_100.mod.dist -f ${read} -i > pggb_100_${read.simpleName}_giraffe.gam
	"""

}

process vgCall {

	publishDir './vg2/', mode: 'copy', overwrite: false

	input:
	file gam from gamfile

	output:
	file '*.vcf'
	file '*.packtable' into packtable
	//path "pggb_100_${read.simpleName}.fastq"

	script:

	"""
	echo "vg pack"
	/home/programs/vg/vg pack -t 6 -x /home/SentericaVariantPipeline/100lines/pggb_100.mod.xg -g ${gam} -o ${gam.simpleName}.pack -d > ${gam.simpleName}.packtable
	/home/programs/vg/vg call -t 6 -d 1 -r /home/SentericaVariantPipeline/100lines/pggb_100.mod.snarls -k ${gam.simpleName}.pack -a -s ${gam.simpleName} /home/SentericaVariantPipeline/100lines/pggb_100.mod.xg > ${gam.simpleName}_allpaths.vcf
	"""
	//echo "vg call against reference"
	// /home/programs/vg/vg call -t 6 -p $refname -d 1 -r ${idxbase}.snarls -k ${out_base}.pack -a -s ${out_base} $xg > ${out_base}.vcf

}

process segcov {

	publishDir './vg2/', mode: 'move', overwrite: false

	input:
	file pack from packtable

	output:
	file '*.cov'

	script:

	"""
	perl /home/SentericaVariantPipeline/pack_table_to_seg_cov.pl -p ${pack} > ${pack.simpleName}.seg.cov
	"""

}

/*
process varCall {

	publishDir './vg2/', mode: 'move', overwrite: false

	input:
	file read from reads_unzippedPGRGvarcall

	output:
	file '*'
	//path "pggb_100_${read.simpleName}.fastq"

	script:
	"""
	/home/SentericaVariantPipeline/vg_giraffe_pipeline.sh pggb_100 ${read} > ${read.simpleName}_log.txt 2>&1
vg giraffe -t 6 -M 20 -g ${idxbase}.gg -H ${idxbase}.gbwt -m ${idxbase}.min -d ${idxbase}.dist -f ${fq1} -i > $gam
	"""

}

*/

/*
sequences = Channel.fromPath('*.fa')
methods = ['regular', 'expresso']
libraries = [ file('PQ001.lib'), file('PQ002.lib'), file('PQ003.lib') ]

process alignSequences {
  input:
  file seq from sequences
  each mode from methods
  each file(lib) from libraries

  """
  t_coffee -in $seq -mode $mode -lib $lib > result
  """
}
*/
