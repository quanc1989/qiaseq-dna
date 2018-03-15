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
RUN mkdir -p /srv/qgen/code && \
    mkdir -p /srv/qgen/bin/downloads && \
    mkdir -p /srv/qgen/data/genome && \
    mkdir -p /srv/qgen/data/annotation && \
    mkdir -p /srv/qgen/example/


################ Update package repository and install dependencies using apt-get ################
RUN apt-get -y update && \
    apt-get -y install r-base

################ Install various version specific 3rd party tools ################
RUN conda install bedtools=2.25.0 htslib=1.3.1 cutadapt=1.10 picard=1.97 snpeff=4.2 bwa=0.7.15
RUN wget https://storage.googleapis.com/qiaseq-dna/lib/ssw.tar.gz \
    https://storage.googleapis.com/qiaseq-dna/lib/fgbio-0.1.4-SNAPSHOT.jar -P /srv/qgen/bin/
RUN cd /srv/qgen/bin/ && \
    tar -xvf ssw.tar.gz    

################ Install python modules ################
## Install some modules with conda
RUN conda install scipy MySQL-python openpyxl pysam=0.9.0
RUN pip install statistics
## Download and install 3rd party libraries
RUN wget https://storage.googleapis.com/qiaseq-dna/lib/py-editdist-0.3.tar.gz https://storage.googleapis.com/qiaseq-dna/lib/sendgrid-v2.2.1.tar.gz -P /srv/qgen/bin/downloads/
RUN cd /srv/qgen/bin/downloads/ && \
    tar -xvf py-editdist-0.3.tar.gz && \
    cd py-editdist-0.3 && \
    /opt/conda/bin/python setup.py install && \
    cd /srv/qgen/bin/downloads/ && \
    tar -xvf sendgrid-v2.2.1.tar.gz && \
    cd sendgrid-python-2.2.1 && \
    /opt/conda/bin/python setup.py install
    
################ R packages ################
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('plyr')"

################ Update openjdk ################
## note : picard gets updated to match jdk version
RUN conda install -c cyclus java-jdk=8.45.14

################ Add latest samtools version for sort by Tag feature ################
RUN wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 -O /srv/qgen/bin/downloads/samtools-1.5.tar.bz2 && \
    cd /srv/qgen/bin/downloads/ && \
    tar -xvf samtools-1.5.tar.bz2 && \
    cd samtools-1.5  && \
    mkdir -p /srv/qgen/bin/samtools-1.5 && \
    ./configure --prefix /srv/qgen/bin/samtools-1.5 && \
    make && \
    make install 

################ Add data directory ################
## Download genome files
RUN wget https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.dict \
         https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.fa.gz -P /srv/qgen/data/genome/
RUN cd /srv/qgen/data/genome && \
    gunzip ucsc.hg19.fa.gz  && \
    ## Index the fasta using samtools
    /srv/qgen/bin/samtools-1.5/bin/samtools faidx /srv/qgen/data/genome/ucsc.hg19.fa && \ 
    ## Run bwa to generate index files 
    /opt/conda/bin/bwa index /srv/qgen/data/genome/ucsc.hg19.fa
    
## Download Annotation files
RUN wget https://storage.googleapis.com/qiaseq-dna/data/annotation/clinvar_20160531.vcf.gz \
         https://storage.googleapis.com/qiaseq-dna/data/annotation/clinvar_20160531.vcf.gz.tbi \
         https://storage.googleapis.com/qiaseq-dna/data/annotation/common_all_20160601.vcf.gz \
    	 https://storage.googleapis.com/qiaseq-dna/data/annotation/common_all_20160601.vcf.gz.tbi \
	 https://storage.googleapis.com/qiaseq-dna/data/annotation/CosmicAllMuts_v69_20140602.vcf.gz \
	 https://storage.googleapis.com/qiaseq-dna/data/annotation/CosmicAllMuts_v69_20140602.vcf.gz.tbi \
	 https://storage.googleapis.com/qiaseq-dna/data/annotation/simpleRepeat_TRF.bed \
	 https://storage.googleapis.com/qiaseq-dna/data/annotation/SR_LC_SL_RepeatMasker.bed \
	 https://storage.googleapis.com/qiaseq-dna/data/annotation/bkg.error.v2.RData \
	 https://storage.googleapis.com/qiaseq-dna/data/annotation/SR_LC_SL.full.bed \
	 https://storage.googleapis.com/qiaseq-dna/data/annotation/simpleRepeat.full.bed \
	  -P /srv/qgen/data/annotation/

## Download annotation using SnpEff command
RUN wget http://downloads.sourceforge.net/project/snpeff/databases/v4_2/snpEff_v4_2_GRCh37.75.zip -P /opt/conda/share/snpeff-4.2-0/
RUN rm -rf /opt/conda/share/snpeff-4.2-0/data/
RUN cd /opt/conda/share/snpeff-4.2-0/ && \
    unzip snpEff_v4_2_GRCh37.75.zip

## The command below is not working anymore because of some certificate issue (debug later)
#RUN /opt/conda/jre/bin/java -jar /opt/conda/share/snpeff-4.2-0/snpEff.jar download GRCh37.75


################ Modules for CNV Analysis ################
## Perl
RUN cpan DateTime
RUN cpan DBI
RUN cpan DBD::SQLite
RUN cpan Env::Path
RUN cpan File::chdir
RUN cpan Getopt::Long::Descriptive
RUN cpan Sort:Naturally
RUN cpan Config::IniFiles
RUN cpan Data::Dump::Color
RUN cpan Data::Table::Excel
RUN cpan Hash::Merge
RUN cpan File::Slurp
## R
RUN Rscript -e "install.packages('MASS')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('gridExtra')"
RUN Rscript -e "install.packages('naturalsort')"
RUN Rscript -e "install.packages('scales')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('extrafont')"
## Annotation file
RUN wget https://storage.googleapis.com/qiaseq-dna/data/annotation/refGene.txt \
         -P /srv/qgen/data/annotation/

################ TVC binaries ################
RUN mkdir -p /srv/qgen/bin/TorrentSuite/
RUN wget https://storage.googleapis.com/qiaseq-dna/lib/TorrentSuite/tmap \
    	 https://storage.googleapis.com/qiaseq-dna/lib/TorrentSuite/tvc \
	 -P /srv/qgen/bin/TorrentSuite/
RUN chmod 775 /srv/qgen/bin/TorrentSuite/tmap /srv/qgen/bin/TorrentSuite/tvc



## Add example fastqs and files
RUN wget https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R1_001.fastq.gz \
    	 https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R2_001.fastq.gz \
	 https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.primers.txt  \
	 https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.roi.bed \
    	 -P /srv/qgen/example/

## Add test files for smCounterv2
RUN mkdir -p /srv/qgen/test_smcounter-v2/
RUN wget https://storage.googleapis.com/qiaseq-dna/test_files/high.confidence.variants.bed \
         https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam \
	 https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.VariantList.long.txt \
	 https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam.bai \
	 -P /srv/qgen/test_smcounter-v2/

################ Update Environment Variables ################
ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python2.7/site-packages/:/srv/qgen/code/qiaseq-dna/
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/srv/qgen/bin/ssw/src/

