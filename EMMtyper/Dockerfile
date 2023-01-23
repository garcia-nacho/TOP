
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

RUN conda install -c conda-forge -c bioconda -c defaults emmtyper
RUN conda install -c conda-forge -c bioconda -c defaults snippy
RUN conda install -c conda-forge ncurses
RUN conda install -c conda-forge -c bioconda -c defaults gubbins
RUN conda install -c conda-forge -c r -c bioconda -c defaults snp-sites
RUN conda install -c bioconda fasttree


USER root  
RUN mkdir -p /Data 

RUN chmod 777 /Data 
USER docker
WORKDIR /Data