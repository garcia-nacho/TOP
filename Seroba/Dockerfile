# base image
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
    && conda create -n seroba 
# install seroba and create database (/seroba-1.0.2/database/)

RUN cd /home/docker && git clone https://github.com/sanger-pathogens/seroba 

RUN /bin/bash -c ". activate seroba && \
    conda install -c bioconda -c conda-forge -c defaults Biopython=1.74 && \
    conda install -c bioconda -c conda-forge -c defaults tbb=2020.3 && \
    conda install -c bioconda -c conda-forge -c defaults kmc=3.1.2rc1 && \
    conda install -c bioconda -c conda-forge -c defaults seroba=1.0.2"
    
RUN /bin/bash -c ". activate seroba && cd /home/docker/seroba && seroba createDBs database/ 71"
USER root
RUN mkdir -p /Data /home/docker/CommonFiles
#COPY CommonFiles/ /home/docker/CommonFiles/
RUN chmod 777 /Data 
USER docker
WORKDIR /Data