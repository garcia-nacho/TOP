# base image
FROM garcianacho/fhibase:v1

LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER docker
WORKDIR /home/docker

# --------------------------------------------------------------------
# Install Miniforge instead of Miniconda latest.
# We need a controlled Python version, not Python 3.13.
# --------------------------------------------------------------------
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p /home/docker/miniforge \
    && rm -f /tmp/miniforge.sh

ENV PATH="/home/docker/miniforge/bin:${PATH}"

# --------------------------------------------------------------------
# Python environment.
# Do NOT install KMC with conda.
# KMC will be installed manually as old SeroBA-style KMC 3.0.0.
# --------------------------------------------------------------------
RUN conda config --set channel_priority strict \
    && mamba create -y -n seroba_env \
        -c conda-forge \
        -c bioconda \
        python=3.7 \
        pip \
        "setuptools<58" \
        "wheel" \
        "cython<3" \
        ariba \
        biopython \
        pyyaml \
        beautifulsoup4 \
        requests \
        pyfastaq \
        nose \
    && conda clean -afy

ENV CONDA_DEFAULT_ENV=seroba_env
ENV PATH="/home/docker/miniforge/envs/seroba_env/bin:/home/docker/miniforge/bin:${PATH}"

# --------------------------------------------------------------------
# Clone GlobalPneumoSeq SeroBA.
# --------------------------------------------------------------------
RUN git clone https://github.com/GlobalPneumoSeq/seroba.git /home/docker/seroba

# --------------------------------------------------------------------
# Install native dependencies manually.
#
# Important:
# - KMC 3.0.0, not 3.2.4
# - MUMmer 3.23
# - Bowtie2 2.3.1 legacy
# - CD-HIT 4.6.8
#
# These match the old SeroBA dependency style.
# --------------------------------------------------------------------
RUN mkdir -p /home/docker/seroba/build/bin

# KMC 3.0.0
RUN mkdir -p /home/docker/seroba/build/kmc_extract \
    && wget -q https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.linux.tar.gz -O /tmp/KMC3.linux.tar.gz \
    && tar -xzf /tmp/KMC3.linux.tar.gz -C /home/docker/seroba/build/kmc_extract \
    && find /home/docker/seroba/build/kmc_extract -type f \( -name "kmc" -o -name "kmc_tools" -o -name "kmc_dump" \) -exec cp {} /home/docker/seroba/build/bin/ \; \
    && chmod +x /home/docker/seroba/build/bin/kmc \
                /home/docker/seroba/build/bin/kmc_tools \
                /home/docker/seroba/build/bin/kmc_dump \
    && rm -f /tmp/KMC3.linux.tar.gz

# MUMmer 3.23
RUN wget -q https://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz -O /tmp/MUMmer3.23.tar.gz \
    && tar -xzf /tmp/MUMmer3.23.tar.gz -C /home/docker/seroba/build \
    && cd /home/docker/seroba/build/MUMmer3.23 \
    && make \
    && rm -f /tmp/MUMmer3.23.tar.gz

# Bowtie2 2.3.1 legacy
RUN wget -q https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.1/bowtie2-2.3.1-legacy-linux-x86_64.zip -O /tmp/bowtie2-2.3.1-legacy.zip \
    && unzip -q -n /tmp/bowtie2-2.3.1-legacy.zip -d /home/docker/seroba/build \
    && rm -f /tmp/bowtie2-2.3.1-legacy.zip

# CD-HIT 4.6.8
RUN wget -q https://github.com/weizhongli/cdhit/archive/V4.6.8.tar.gz -O /tmp/cdhit-4.6.8.tar.gz \
    && tar -xzf /tmp/cdhit-4.6.8.tar.gz -C /home/docker/seroba/build \
    && cd /home/docker/seroba/build/cdhit-4.6.8 \
    && make \
    && rm -f /tmp/cdhit-4.6.8.tar.gz

# --------------------------------------------------------------------
# Native tools first in PATH.
# This prevents conda-installed versions from shadowing the old tools.
# --------------------------------------------------------------------
ENV PATH="/home/docker/seroba/build/bin:/home/docker/seroba/build/MUMmer3.23:/home/docker/seroba/build/bowtie2-2.3.1-legacy:/home/docker/seroba/build/cdhit-4.6.8:/home/docker/miniforge/envs/seroba_env/bin:/home/docker/miniforge/bin:${PATH}"

# --------------------------------------------------------------------
# Install Python packages that SeroBA expects.
# Keep pymummer old.
# Use --no-build-isolation to avoid modern build-tool surprises.
# --------------------------------------------------------------------
RUN python -m pip install --no-build-isolation "pymummer==0.10.3"

# Install SeroBA itself, but do not let pip reinstall dependencies.
RUN cd /home/docker/seroba \
    && python -m pip install --no-deps .

# --------------------------------------------------------------------
# Create SeroBA database.
# --------------------------------------------------------------------
RUN cd /home/docker/seroba \
    && seroba createDBs /home/docker/seroba/database/ 71

# --------------------------------------------------------------------
# Sanity checks.
# These should show /home/docker/seroba/build/bin/kmc.
# --------------------------------------------------------------------
RUN which seroba \
    && which kmc \
    && which kmc_tools \
    && which nucmer \
    && which bowtie2 \
    && which cd-hit-est \
    && echo "CD-HIT path: $(which cd-hit-est)" \
    && seroba version \
    && seroba --help

# --------------------------------------------------------------------
# Critical KMC sanity test.
# If this fails, the image should fail during build.
# --------------------------------------------------------------------
RUN cd /tmp \
    && rm -rf tiny.fq tiny_kmc* kmc_tmp \
    && mkdir -p kmc_tmp \
    && printf '@r1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' > tiny.fq \
    && kmc -k21 -ci1 -m2 -t1 tiny.fq tiny_kmc kmc_tmp \
    && ls -lh tiny_kmc.kmc_pre tiny_kmc.kmc_suf

# Optional. Uncomment only if you want build-time tests.
# RUN cd /home/docker/seroba && python setup.py test

USER root
RUN mkdir -p /Data /home/docker/CommonFiles \
    && chmod 777 /Data

USER docker
WORKDIR /Data
