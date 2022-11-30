# base image
FROM garcianacho/fhibase:v1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"
ARG PROKKA_VER="1.14.5"

#Rflavoured from staph-b docker implementation of prokka

# install dependencies
RUN apt-get update && apt-get -y --no-install-recommends install \
 bzip2 \
 gzip \
 wget \
 perl \
 less \
 libdatetime-perl \
 libxml-simple-perl \
 libdigest-md5-perl \
 default-jre \
 bioperl \
 hmmer \
 zlib1g-dev \
 python \
 liblzma-dev \
 libbz2-dev \
 xz-utils \ 
 curl \
 g++ \
 cpanminus \
 make && apt-get autoclean && rm -rf /var/lib/apt/lists/*

ENV PATH="$PATH:/bedtools2/bin"


# install bedtools 2.29.0 since >=2.27.0 is required for barrnap and the apt-get package is 2.25.0
# dependencies required for bedtools: zlib1g-dev python liblzma-dev libbz2-dev xz-utils curl g++
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz && \
  tar -zxf bedtools-2.29.0.tar.gz && \
  rm bedtools-2.29.0.tar.gz && \
  cd bedtools2 && \
  make


# install barrnap
RUN wget https://github.com/tseemann/barrnap/archive/0.9.tar.gz && \
    tar -zxf 0.9.tar.gz && \
    rm 0.9.tar.gz && \
    cd barrnap-0.9 && \
    make test

# to fix error when running Prokka 1.14.5: '"uniq" is not exported by the List::Util module'
# error introduced while updating dockerfile from prokka 1.14.0 to 1.14.5
RUN cpanm List::Util

# download prokka and make /data
RUN wget https://github.com/tseemann/prokka/archive/v${PROKKA_VER}.tar.gz && \
    tar -xzf v${PROKKA_VER}.tar.gz && \
    rm -rf v${PROKKA_VER}.tar.gz && \
    mkdir /data

# install tbl2asn manually since the one included with prokka is expired. 
# Probably will have to do again in Dec 2021 unless Torsten removes it from prokka
RUN wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz -O linux64.tbl2asn.gz && \
    gunzip linux64.tbl2asn.gz && \
    mv linux64.tbl2asn /usr/bin/tbl2asn

RUN chmod 755 /usr/bin/tbl2asn && mkdir -p /Data /home/docker/CommonFiles && chmod 777 /Data 

# set PATH to prioritize /usr/bin/tbl2asn (and other paths) before /prokka-[VERSION]/binaries*/tbl2asn
# second line added so that parallel is forced into PATH (wasn't working in Singularity without explicitly setting this)
ENV PATH="$PATH:/prokka-${PROKKA_VER}/bin:\
/prokka-${PROKKA_VER}/binaries/common:\
/barrnap-0.9/bin"

# set perl locale settings so barrnap works in singularity image
ENV LC_ALL=C

# index dbs and list what dbs were setup; set working dir to /data
RUN prokka --setupdb && prokka --listdb
USER docker
WORKDIR /Data
