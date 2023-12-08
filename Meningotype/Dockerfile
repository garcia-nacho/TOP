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

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends libmoo-perl liblist-moreutils-perl libjson-perl curl


USER docker

RUN conda install -c bioconda any2fasta \
    && conda install -c bioconda mamba \
    && conda install -c conda-forge biopython \
    && conda install -c bioconda ispcr \ 
    && conda install -c bioconda blast
RUN conda install -c bioconda mlst

RUN pip install --user git+https://github.com/MDU-PHL/meningotype.git

USER root
RUN ln -s /home/docker/.local/bin/* /usr/bin
RUN apt-get update \
    && apt-get install -y --no-install-recommends file
USER docker