#!/bin/bash

#Merging from WW
#Nextflow
#Cleaning and organizing
#docker run -it --rm -v $(pwd):/Data -v /media/nacho/Data/kraken2_standard_20220926/:/Kraken2DB garcianacho/top bash

conda create -n top_nf -y -c bioconda nextflow

#If ${1} doesnt exist go to pwd
echo "Downloading Kraken database"
echo ""
wget -O ${1}/krakenDB https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz

echo "Unpacking Kraken database"
echo ""
tar -xvzf ${1}/krakenDB https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz

echo "Installing docker images"
echo ""

docker get garcianacho/top:spades
docker get garcianacho/top:seroba
docker get garcianacho/top:hicap
docker get garcianacho/top:abricate
docker get garcianacho/top:prokka

docker build -t garcianacho/top:spades Spades/
docker build -t garcianacho/top:seroba Seroba/
docker build -t garcianacho/top:hicap HiCap/
docker build -t garcianacho/top:abricate Abricate/
docker build -t garcianacho/top:prokka Prokka/