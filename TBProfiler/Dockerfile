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
    && conda config --add channels conda-forge

RUN conda install -c conda-forge mamba
RUN mamba create -y -n tb-profiler tb-profiler
USER root
RUN Rscript -e "install.packages(c('readxl','writexl','jsonlite'))"
RUN /bin/bash -c ". activate tb-profiler && tb-profiler update_tbdb --branch tbdb"
COPY Code/ /home/docker/Code/
USER docker
WORKDIR /Data
