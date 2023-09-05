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
RUN pip3 install pandas pyfaidx biopython pyahocorasick \
    && cd /home/docker && git clone https://github.com/leosanbu/pyngSTar 

COPY DB/ /home/docker/pyngSTar/pyngSTarDB_050923
COPY Code /home/docker/Code
USER root
RUN chmod -R +rwx /home/docker/pyngSTar/pyngSTarDB_050923/* \
    && chmod 777 /home/docker/pyngSTar/pyngSTarDB_050923

RUN chmod -R +rwx /home/docker/Code/* \
    && ln /home/docker/Code/ngstar.sh /usr/bin/ngstar.sh
USER docker
RUN cd /home/docker/pyngSTar/ && python3 pyngSTar.py -p pyngSTarDB_050923 -u
CMD /usr/bin/ngstar.sh