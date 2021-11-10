params.index = "/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/3SentericaSRAAccessions.txt"

params.project = "SRR12660755"

params.resultdir = '/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/results'

params.outputdir = "/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/work"

params.qualevaldir = "/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/FastQC"

LT2ref = Channel.fromPath( '/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/AE006468.2.fasta' )

projectSRId = params.project

int threads = Runtime.getRuntime().availableProcessors()

process getSRA {
		
	publishDir params.resultdir, mode: 'copy'

	cpus 1

	input:
	path 'accIDs' from params.index
	val projectID from projectSRId
	
	output:
	file 'sra.txt' into sraIDs
	
	script:
	"""
	tail accIDs -n +2 > sra.txt
	"""
}

sraIDs.splitText().map { it -> it.trim() }.set { singleSRAId }
singleSRAId.into { singleSRAIdFQD; singleSRAIdMR }

process fastqDump {

	publishDir params.resultdir, mode: 'copy'

	input:
	val id from singleSRAIdFQD

	output:
	file '*.fastq.gz' into reads

	//Note that the number of spot IDs should be determined based on read length for a dataset in order to achieve 60X genome coverage

	script:
	"""
	fastq-dump --clip --skip-technical --readids --maxSpotId 2000000 --gzip ${id}
	"""

	// --read-filter pass
	// --split-files
	//parallel-fastq-dump --sra-id $id --threads ${task.cpus} --gzip
}

reads.into { reads_dashing; reads_fastqc; reads_gunzip; }

process dashingSketch {

	publishDir '/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/dashing', mode: 'move', overwrite: false

	cpus threads

	input:
	file read from reads_dashing
	
	output:
	file '*.hll'
	
	script:
	readName = read.toString() - ~/(\.fastq\.gz)?$/
	
	"""
	/home/chrisgaby/program_executables/dashing/dashing sketch --nthreads $threads -S 10 -k 16 $read
	"""
}

process fastQC {

	//publishDir params.resultdir, mode: 'move'

	cpus threads

	input:
	file read from reads_fastqc
	
	//output:
	//file '*.zip' into publishDir
	
	script:
	
	"""
	~/program_executables/fastqc_v0.11.8/FastQC/fastqc -t $threads -f fastq --extract -o $params.qualevaldir $read
	"""
}

process indexRef {

	publishDir '/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/', mode: 'move', overwrite: false

	input:
	file ref from LT2ref

	output:
	file '*'
	//into bwa_index

	script:

	"""
	bwa index $ref
	"""
}

process gunzipReads {

	publishDir params.resultdir, mode: 'copy', overwrite: false

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

process mapReads {

	publishDir '/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/mappedToRef/', mode: 'copy', overwrite: false

	cpus threads

	input:
	file read from reads_unzipped
	//file '*' from bwa_index.first()
	//val id from singleSRAIdMR

	output:
	file '*.bam' into bamFile

	script:

	"""
	bwa mem /home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/AE006468.2.fasta /home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/results/${read} -t $threads | samtools sort -o ${read.simpleName}ChromMappedToLT2ref.bam
	samtools index ${read.simpleName}ChromMappedToLT2ref.bam
	"""
}

bamFile.into { bamFile_stats; bamFile_mpileup; bamFile_freebayes }

process mapStats {

	publishDir '/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/mapStats/', mode: 'move', overwrite: false

	input:
	file bamFile from bamFile_stats

	output:
	file "*.txt"

	script:

	"""
	samtools flagstat ${bamFile} > ${bamFile.simpleName}Stats.txt
	"""

}

process mpileup {

	//publishDir '/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/mpileup', mode: 'move', overwrite: false

	input:
	file bamFile from bamFile_mpileup

	//output:
	//file "*.vcf"

	script:

	"""
	bcftools mpileup -Ou -f /home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/AE006468.2.fasta ${bamFile} | bcftools call -mv -Ov -o /home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/mpileup/${bamFile.simpleName}_MP.vcf
	"""
	//samtools mpileup -f /home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/AE006468.2.fasta -g ${bamFile} | bcftools call -mv - > /home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/mpileup/${bamFile.simpleName}_MP.vcf
}

process freebayes {

	//publishDir '/home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/freebayes', mode: 'move', overwrite: false

	input:
	file bamFile from bamFile_freebayes

	//output:
	//file "*.vcf"

	script:
	"""
	freebayes -f /home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/AE006468.2.fasta -p 1 ${bamFile} > /home/chrisgaby/Documents/USDA/GoogleCloudSimpleNextflow/freebayes/${bamFile.simpleName}_FB.vcf
	"""

}

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
