## Docker Image with Variant and SNP Calling Pipelines

An intial docker image was assembled by John Christian Gaby on June 28, 2021. Subsequently, Brian Nadon and Justin Vaughn modified the variant calling pipeline to utilize the Giraffe aligner of vg tools, and John Chrisitan Gaby updated the Dockerfile to include dependencies required for the updated variant calling pipeline as well as additional processes that were incorporated into the Nextflow pipeline script, such as FastQC evaluation of SRA read datasets and fastq-dump (from SRA tools) download of read datasets from the SRA.

This folder, named 'dockerbuild', contains the Dockerfile needed to build a Docker image for variant and SNP calling of *S. enterica* genome sequence data.

The file named 'Dockerfile' is a series of commands that Docker will use to assemble the files into a working folder as well as install programs and download files from the Internet into the Docker image.

The working folder is '/home/SentericaVariantPipeline/'

The working folder in the image will contain the reference graph and the *S. enterica* LT2 genome, which are needed for variant and SNP calling, respectively. 

Also in the working folder is a script for running a demo of the variant calling pipeline:

'Graph_align_pipeline_test.sh'

And there is a second script for running a demo of the SNP calling pipeline:

'SNPcall.sh'

The output of the variant calling pipeline is a '.segcov' file, whereas the output of the SNP calling pipeline is a '.vcf' file.

Installed executables are located as follows:

/root/miniconda3/bin/conda
/root/miniconda3/bin/bwa
/root/miniconda3/bin/GraphAligner
/home/programs/samtools_1.12/bin/samtools
/home/programs/bcftools_1.12/bin/bcftools
/home/programs/vg/vg

To build the Docker image, first change directories to the 'dockerbuild' directory with the 'dockerbuild' file, and then run the following command that will run the 'dockerbuild' file to construct the image:

'docker build -t chrisgaby/vcall:latest -t chrisgaby/vcall:0.1 .'

Provided that Docker succeeds in constructing the image, then one should be able to list the images with the following command:

'docker images'

And there should be an image called 'chrisgaby/vcall', with both the 'latest' and '0.1' versions having an identical image ID. If the image exists, then one should be able to run an interactive container of the image with the following command:

'docker run -it chrisgaby/vcall bash'

This should switch to the container terminal, and if one types the bash command 'pwd' to print the working directory, then one should be in the folder '/home/SentericaVariantPipeline/'.

Then, one may run the demo script for the variant calling pipeline by typing the following into the terminal:

'./Graph_align_pipeline_test.sh'

it will take ~10 minutes to run. Then, the SNP calling demo pipeline may be run by typing the following into the terminal:

'./SNPcall.sh'

While the demo scripts process just one read dataset from the Sequence Read Archive (SRA), the goal is to create a Nextflow pipeline that uses this docker image to process the >300,000 *S. enterica* isolates whose genomes have been sequenced and are available from the SRA.

The docker image is 6.4 GB in size. However, about 4.5 GB of that size is attributable to the 2 reference graph files, 'pggb_250lines.gfa' and 'pggb_250lines.xg', that total 4.5 GB. For the scaled-up run to analyze the >300,000 *S. enterica* SRA datasets, these files for the reference graph may be stored on a single volume in the cloud rather than being distributed with the docker image, thus reducing the docker image size.
