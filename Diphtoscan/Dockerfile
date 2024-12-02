
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

RUN conda install -c bioconda mash
RUN conda install -c bioconda blast
RUN conda install -c bioconda ncbi-amrfinderplus
RUN conda install --solver=classic conda-forge::conda-libmamba-solver conda-forge::libmamba conda-forge::libmambapy conda-forge::libarchive
RUN conda install -c bioconda hmmer
RUN conda install pandas
RUN conda install biopython

USER root  

COPY diphtoscan /home/docker/diphtoscan
RUN /home/docker/diphtoscan/script/update_tools.sh
RUN chmod -R 777 /home/docker/diphtoscan

USER docker
RUN amrfinder_update -d /home/docker/diphtoscan/data/resistance/

WORKDIR /Data
CMD ["sh", "-c", "/home/docker/diphtoscan/Diphtorunner.sh"]