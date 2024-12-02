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

RUN conda install -y -c conda-forge mamba
RUN mamba install -y -c conda-forge -c bioconda -c defaults snippy
RUN mamba install -y -c conda-forge ncurses
RUN mamba install -y -c conda-forge -c bioconda -c defaults gubbins
RUN mamba install -y -c conda-forge -c r -c bioconda -c defaults snp-sites
RUN mamba install -y -c bioconda fasttree
RUN mamba install -y -c bioconda bowtie2

USER root
# Install OpenJDK-8
#RUN apt-get update && \
#    apt-get install -y software-properties-common && \
#    add-apt-repository ppa:linuxuprising/java && \
#    apt-get update && \
#    apt-get install -y openjdk-8-jdk && \
#    apt-get install -y ant && \
#    apt-get clean;

RUN apt-get update && \
    apt-get install -y default-jdk

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

RUN Rscript -e "install.packages(c('sp','seqinr', 'ggplot2', 'writexl', 'ape'), repos = 'http://cran.us.r-project.org')"

USER docker

COPY Refs /home/docker/references
COPY snpEff /home/docker/snpEff
COPY code /home/docker/code
COPY bt2 /home/docker/bt2
COPY blastdb /home/docker/blastdb
COPY profiles /home/docker/profiles

WORKDIR /Data
CMD ["bash", "-c", "/home/docker/code/bpe_mlst.sh"]