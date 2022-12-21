
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
    && conda create -n hicap

RUN /bin/bash -c ". activate hicap && \
   conda install -c bioconda -c conda-forge hicap"  
USER root  
RUN mkdir -p /Data /home/docker/Code
COPY Code/ /home/docker/Code
RUN chmod -R +rwx /home/docker/Code/* \
    && chmod 777 /Data 
USER docker
WORKDIR /Data