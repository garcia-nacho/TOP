FROM ghcr.io/garcia-nacho/top_spades:v1.1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

USER root

RUN mkdir -p /Data /home/docker/CommonFiles
COPY CommonFiles/ /home/docker/CommonFiles/
RUN chmod -R +rwx /home/docker/CommonFiles/* \
    && chmod 777 /Data 

RUN cd /home/docker/CommonFiles/ && tar -xf blast_mlst.tar.gz

USER docker
WORKDIR /Data