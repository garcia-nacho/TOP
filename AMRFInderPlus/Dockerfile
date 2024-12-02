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

RUN conda install -y blast hmmer libcurl
RUN cd /home/docker/ && mkdir amrfinder && cd amrfinder
RUN cd /home/docker/amrfinder && wget https://github.com/ncbi/amr/releases/download/amrfinder_v3.12.8/amrfinder_binaries_v3.12.8.tar.gz && tar -xf amrfinder_binaries_v3.12.8.tar.gz
RUN /home/docker/amrfinder/amrfinder -u

USER docker
WORKDIR /Data