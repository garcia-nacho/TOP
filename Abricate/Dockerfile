#Fully-loaded version derived from StatPH-B Fl but based on fhibase
FROM garcianacho/fhibase:v1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"


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
    && conda install -y -c conda-forge mamba

RUN mamba install -c conda-forge -c bioconda -c defaults abricate
USER root
RUN mkdir /Data
RUN chmod 777 /Data 

USER root
#RUN abricate-get_db --db ncbi --force &&\
#    abricate-get_db --db vfdb --force &&\
#    abricate-get_db --db plasmidfinder --force
COPY Code/ /home/docker/Code
RUN chmod -R +rwx /home/docker/Code/*

COPY vfdb/ /home/docker/miniconda3/db/vfdb2
RUN chmod 777 /home/docker/miniconda3/db/vfdb2/sequences.gz
RUN cd /home/docker/miniconda3/db/vfdb2/ && gzip -d sequences.gz
RUN abricate --setupdb
WORKDIR /Data
