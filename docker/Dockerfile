FROM ubuntu:xenial

# create working directories
RUN mkdir /home/SentericaVariantPipeline/
RUN mkdir /home/programs/

# set the working directory
WORKDIR /home/SentericaVariantPipeline/

# install variant calling files and programs
# copy the reference graph files into the image
COPY pggb_250lines.gfa pggb_250lines.xg ./

# update the aptitude package manager repositories
RUN apt-get update

# install needed software packages wget and python3
RUN apt-get -y install wget python3

# download a required script
RUN wget https://github.com/brianabernathy/gfa_var_genotyper/raw/b01900f2fac6370da3df6070d82fecdbfd3e3c80/pack_table_to_seg_cov.pl

# install vg
RUN mkdir /home/programs/vg/ && cd /home/programs/vg/ && wget https://github.com/vgteam/vg/releases/download/v1.30.0/vg && chmod +x vg

# install miniconda and use it to install graphaligner and bwa
RUN cd /home/programs/ && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash ./Miniconda3-latest-Linux-x86_64.sh -b && rm ./Miniconda3-latest-Linux-x86_64.sh && /root/miniconda3/bin/conda install -c bioconda graphaligner=1.0.11 && /root/miniconda3/bin/conda install -c bioconda bwa=0.7.17

# install SNP calling files and programs
# copy fastq reads from the SRA dataset for Salmonella enterica subsp. enterica serovar Tennessee
RUN wget https://sra-download.ncbi.nlm.nih.gov/traces/sra79/SRZ/014464/SRR14464606/87712.fastq -O SRR14464606_1.fastq && wget https://sra-download.ncbi.nlm.nih.gov/traces/sra79/SRZ/014464/SRR14464606/87713.fastq -O SRR14464606_2.fastq

# copy the S. enterica LT2 reference genome and its plasmid
COPY AE006468.2.fasta AE006471.2.fasta ./

# install dependency to unpack bz2 files prior to its use in samtools install
RUN apt-get install libbz2-dev lbzip2 -y

# download and unpack samtools
RUN mkdir /home/programs/samtools/ && cd /home/programs/samtools/ && wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && tar -xf *.bz2 && rm *.bz2

# install samtools dependencies
RUN apt-get install gcc libncurses5-dev zlib1g-dev liblzma-dev libcurl4-openssl-dev libssl-dev build-essential -y

# install samtools
RUN cd /home/programs/samtools/samtools-1.12/ && mkdir /home/programs/samtools_1.12/ && ./configure --prefix=/home/programs/samtools_1.12/ && make && make install && rm /home/programs/samtools/ -r

# install bcftools
RUN mkdir /home/programs/bcftools/ && cd /home/programs/bcftools/ && wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 && tar -xf *.bz2 && rm *.bz2 && cd /home/programs/bcftools/bcftools-1.12/ && mkdir /home/programs/bcftools_1.12/ && ./configure --prefix=/home/programs/bcftools_1.12/ && make && make install && rm /home/programs/bcftools -r

# copy the original and modified version of the variant calling pipeline script
COPY Graph_align_pipeline.sh Graph_align_pipeline_test.sh ./

# make the demo variant calling shell script executable
RUN chmod +x Graph_align_pipeline_test.sh

# copy the SNP calling script
COPY SNPcall.sh ./

# make the demo SNP calling shell script executable
RUN chmod +x SNPcall.sh
