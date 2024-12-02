
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
RUN conda create -n emmtip
RUN /bin/bash -c ". activate emmtip && \
    pip install git+https://github.com/admiralenola/emm_typing"

USER root 
RUN Rscript -e "install.packages(c( 'seqinr'),repos = 'http://cran.us.r-project.org')" 
RUN mkdir -p /Data 
RUN mkdir -p /emmdb07022024 
COPY emmdb.nhr /emmdb07022024/emmdb.tfa.nhr
COPY emmdb.nin /emmdb07022024/emmdb.tfa.nin
COPY emmdb.nsq /emmdb07022024/emmdb.tfa.nsq
COPY alltrimmed.fa /emmdb07022024/emmdb.tfa
COPY EMM_Extraction.R /home/docker/EMM_Extraction.R
RUN chmod -R +rx /emmdb07022024/*

RUN chmod 777 /Data 
RUN chmod +x /home/docker/EMM_Extraction.R

USER docker
WORKDIR /Data
CMD ["bash", "-c", "Rscript /home/docker/EMM_Extraction.R"]