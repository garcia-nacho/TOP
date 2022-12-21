#Fully-loaded version derived from StatPH-B Fl but based on fhibase
FROM garcianacho/fhibase:v1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"
ARG ABRICATE_VER="1.0.1"

# For easy upgrade later. ARG varibles only persist during docker image build time
ARG ANY2FASTA_VERSION="0.4.2"

LABEL description="Mass screening of contigs for AMR or virulence genes"
LABEL website="https://github.com/tseemann/abricate"
LABEL license="https://github.com/tseemann/abricate/blob/master/LICENSE"
LABEL maintainer="Nacho Garcia"
LABEL maintainer.email="iggl@fhi.no"

# install dependencies
# removed: emboss
# ncbi-blast+ version in apt for ubuntu:focal = v2.9.0 
RUN apt-get update && apt-get install -y --no-install-recommends \
  bioperl \
  gzip \
  unzip \
  liblist-moreutils-perl \
  libjson-perl \
  libtext-csv-perl \
  libfile-slurp-perl \
  liblwp-protocol-https-perl \
  libwww-perl \
  libpath-tiny-perl \
  git \
  ncbi-blast+ \
  wget && \
  apt-get autoclean && rm -rf /var/lib/apt/lists/*

# get any2fasta
RUN wget https://github.com/tseemann/any2fasta/archive/refs/tags/v${ANY2FASTA_VERSION}.tar.gz && \
  tar -xvf v${ANY2FASTA_VERSION}.tar.gz && \
  rm v${ANY2FASTA_VERSION}.tar.gz && \
  cd any2fasta-${ANY2FASTA_VERSION} && \
  chmod +x any2fasta && \
  cp any2fasta /usr/local/bin

# download abricate; make /data; check dependencies
RUN wget https://github.com/tseemann/abricate/archive/v${ABRICATE_VER}.tar.gz && \
  tar -zxvf v${ABRICATE_VER}.tar.gz && \
  rm -rf v${ABRICATE_VER}.tar.gz && \
  mkdir /Data && \
  /abricate-${ABRICATE_VER}/bin/abricate --check

# set $PATH
# set perl locale settings for singularity compatibility
ENV PATH="/abricate-${ABRICATE_VER}/bin:\
$PATH"\
  LC_ALL=C

# add custom database; rename FASTA file as 'sequences'; index database with --setupdb
ADD AcinetobacterPlasmidTyping.fasta /
RUN mkdir -v /abricate-${ABRICATE_VER}/db/AcinetobacterPlasmidTyping && \
  cp -v /AcinetobacterPlasmidTyping.fasta /abricate-${ABRICATE_VER}/db/AcinetobacterPlasmidTyping/sequences && \
  abricate --setupdb && \
  abricate --list

ADD Hi_Fts1.fasta /
RUN mkdir -v /abricate-${ABRICATE_VER}/db/HinfFtsI && \
  cp -v /Hi_Fts1.fasta /abricate-${ABRICATE_VER}/db/HinfFtsI/sequences && \
  abricate --setupdb && \
  abricate --list

ADD Hi_GyrSubA.fasta /
RUN mkdir -v /abricate-${ABRICATE_VER}/db/HinfGyrSubA && \
  cp -v /Hi_GyrSubA.fasta /abricate-${ABRICATE_VER}/db/HinfGyrSubA/sequences && \
  abricate --setupdb && \
  abricate --list

ADD Hi_TopoIVsubA.fasta /
RUN mkdir -v /abricate-${ABRICATE_VER}/db/HinfTopoIVsubA && \
  cp -v /Hi_TopoIVsubA.fasta /abricate-${ABRICATE_VER}/db/HinfTopoIVsubA/sequences && \
  abricate --setupdb && \
  abricate --list

RUN chmod -R 777 /abricate-${ABRICATE_VER}/
# final working directory is /data
RUN mkdir -p /Data /home/docker/Code
COPY Code/ /home/docker/Code
RUN chmod -R +rwx /home/docker/Code/* \
    && chmod 777 /Data 
USER docker
WORKDIR /Data
