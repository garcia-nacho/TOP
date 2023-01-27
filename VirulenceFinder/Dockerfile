FROM garcianacho/fhibase:v1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

### RUN set -ex; \
USER root 
RUN apt-get update -qq; \
    apt-get install -y -qq git \
    apt-utils \
    wget \
    python3-pip \
    ncbi-blast+ \
    libz-dev \
    ; \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*;
USER docker
# Install python dependencies
RUN pip3 install -U biopython==1.73 tabulate cgecore

# Install kma
RUN cd /home/docker && \
    git clone --branch 1.0.1 --depth 1 https://bitbucket.org/genomicepidemiology/kma.git; \
    cd kma && make

USER root

RUN mv /home/docker/kma/kma* /bin/ && mkdir /Data
COPY virulencefinder.py /usr/src/virulencefinder.py
COPY db /database

RUN chmod 755 /usr/src/virulencefinder.py
RUN chmod -R +rwx /database/* \
    && chmod 777 /Data 
ENV PATH $PATH:/usr/src

USER docker
WORKDIR /Data
