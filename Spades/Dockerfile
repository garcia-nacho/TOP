FROM garcianacho/fhibase:v1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"
USER docker
RUN cd /home/docker \
    && wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /home/docker/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -c bioconda spades \
    && conda install -c bioconda samtools\
    && conda install -c bioconda seqkit \
    && conda update -n base -c defaults conda \
    && conda install -c bioconda bedtools \
    && conda install -c bioconda bowtie2 \
    && conda install -c bioconda varscan \
    && conda install -c bioconda trimmomatic \
    && conda install -c bioconda fastqc \
    && conda install -c bioconda multiqc \
    && conda install nextflow \
    && ln -s /home/docker/miniconda3/lib/libcrypto.so.1.1 /home/docker/miniconda3/lib/libcrypto.so.1.0.0    

USER root

RUN Rscript -e "install.packages(c( 'seqinr', 'ggplot2',  'reshape2', 'ggpubr', 'readxl','tidyverse','writexl', 'devtools', 'data.table','digest'))"
ARG K2VER="2.1.2"
RUN apt-get update && apt-get -y --no-install-recommends install \
 cpanminus \
 ncbi-blast+ \
 libncurses5 && \
 rm -rf /var/lib/apt/lists/* && apt-get autoclean

RUN cpanm Getopt::Std

# DL Kraken2, unpack, and install
RUN wget https://github.com/DerrickWood/kraken2/archive/v${K2VER}.tar.gz && \
 tar -xzf v${K2VER}.tar.gz && \
 rm -rf v${K2VER}.tar.gz && \
 cd kraken2-${K2VER} && \
 ./install_kraken2.sh . && \
 mkdir /kraken2-db

RUN Rscript -e "install.packages(c('jsonlite', 'R.utils'))"

RUN mkdir -p /Data /home/docker/CommonFiles
COPY CommonFiles/ /home/docker/CommonFiles/
RUN chmod -R +rwx /home/docker/CommonFiles/* \
    && chmod 777 /Data 
#Unpack Kraken2 db
USER docker
WORKDIR /Data
