####################################################################################
########        This is a Dockerfile to describe QIAGEN's read processing  #########
########        runtime framework for spe-dna panels                       #########
####################################################################################

# Using a biocontainer base image
# Please see below for further details : 
# https://github.com/BioContainers/containers/blob/master/biocontainers/Dockerfile
FROM biocontainers/biocontainers:latest

MAINTAINER Raghavendra Padmanabhan <raghavendra.padmanabhan@qiagen.com>

################ Create appropriate directory structure for code to run ################
USER root
RUN mkdir -p /srv/qgen/code
RUN mkdir -p /srv/qgen/bin/downloads
RUN mkdir -p /srv/qgen/data/genome
RUN mkdir -p /srv/qgen/data/annotation
RUN mkdir -p /srv/qgen/example

################ Update package repository ################
RUN apt-get -y update

################ Install various version specific 3rd party tools ################
RUN conda install bedtools=2.25.0
RUN conda install htslib=1.3.1
RUN conda install cutadapt=1.10
RUN conda install picard=1.97
RUN conda install snpeff=4.2
RUN conda install bwa=0.7.15
RUN conda install gatk=3.7

################ Install python modules ################
## Install some modules with conda
RUN conda install scipy
RUN conda install MySQL-python
RUN conda install openpyxl
RUN conda install pysam=0.9.0
## Download and install 3rd party libraries
ADD https://storage.googleapis.com/qiaseq-dna/lib/py-editdist-0.3.tar.gz /srv/qgen/bin/downloads/
RUN cd /srv/qgen/bin/downloads/ && \
    tar -xvf py-editdist-0.3.tar.gz && \
    cd py-editdist-0.3 && \
    /opt/conda/bin/python setup.py install

ADD https://storage.googleapis.com/qiaseq-dna/lib/sendgrid-v2.2.1.tar.gz /srv/qgen/bin/downloads/
RUN cd /srv/qgen/bin/downloads/ && \
    tar -xvf sendgrid-v2.2.1.tar.gz && \
    cd sendgrid-python-2.2.1 && \
    /opt/conda/bin/python setup.py install
	    
ADD https://storage.googleapis.com/qiaseq-dna/lib/ssw.tar.gz /srv/qgen/bin/
RUN cd /srv/qgen/bin/ && \
    tar -xvf ssw.tar.gz    
ADD https://storage.googleapis.com/qiaseq-dna/lib/fgbio-0.1.4-SNAPSHOT.jar /srv/qgen/bin/downloads/

################ Install R ################
RUN apt-get -y install r-base

################ Update openjdk ################
## note : picard gets updated to match jdk version
RUN conda install -c cyclus java-jdk=8.45.14

################ Add latest samtools version for sort by Tag feature ################
RUN wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 -O /srv/qgen/downloads/samtools-1.5.tar.bz2 && \
    cd /srv/qgen/downloads/ && \
    tar -xvf samtools-1.5.tar.bz2 && \
    cd samtools-1.5  && \
    mkdir -p /srv/qgen/bin/samtools-1.5 && \
    ./configure --prefix /srv/qgen/bin/samtools-1.5 && \
    make && \
    make install 

################ Add data directory ################
## Download genome files
ADD https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.dict /srv/qgen/data/genome/
ADD https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.fa.gz /srv/qgen/data/genome/
RUN cd /srv/qgen/data/genome && \
    gunzip ucsc.hg19.fa.gz 
## Run bwa to generate index files 
RUN /opt/conda/bin/bwa index /srv/qgen/data/genome/ucsc.hg19.fa

## Download Annotation files
ADD https://storage.googleapis.com/qiaseq-dna/data/annotation/clinvar_20160531.vcf.gz /srv/qgen/data/annotation/
ADD https://storage.googleapis.com/qiaseq-dna/data/annotation/clinvar_20160531.vcf.gz.tbi /srv/qgen/data/annotation/
ADD https://storage.googleapis.com/qiaseq-dna/data/annotation/common_all_20160601.vcf.gz /srv/qgen/data/annotation/
ADD https://storage.googleapis.com/qiaseq-dna/data/annotation/common_all_20160601.vcf.gz.tbi /srv/qgen/data/annotation/
ADD https://storage.googleapis.com/qiaseq-dna/data/annotation/CosmicAllMuts_v69_20140602.vcf.gz /srv/qgen/data/annotation/
ADD https://storage.googleapis.com/qiaseq-dna/data/annotation/CosmicAllMuts_v69_20140602.vcf.gz.tbi /srv/qgen/data/annotation/
ADD https://storage.googleapis.com/qiaseq-dna/data/annotation/simpleRepeat_TRF.bed /srv/qgen/data/annotation/
ADD https://storage.googleapis.com/qiaseq-dna/data/annotation/SR_LC_SL_RepeatMasker.bed /srv/qgen/data/annotation/

## Download annotation using SnpEff command
RUN /opt/conda/jre/bin/java -jar /opt/conda/share/snpeff-4.2-0/snpEff.jar download GRCh37.75

## Add example fastqs and files
RUN wget https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R1_001.fastq.gz -O /srv/qgen/example/NEB_S2_L001_R1_001.fastq.gz
RUN wget https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R2_001.fastq.gz -O /srv/qgen/example/NEB_S2_L001_R2_001.fastq.gz
RUN wget https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.primers.txt -O /srv/qgen/example/DHS-101Z.primers.txt
RUN wget https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.roi.bed -O /srv/qgen/example/DHS-101Z.roi.bed

################ Update Environment Variables ################
ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python2.7/site-packages/:/srv/qgen/code/qiaseq-dna/
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/srv/qgen/bin/ssw/src/
