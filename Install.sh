#!/bin/bash

conda create -n top_nf -y -c bioconda nextflow
source activate top_nf
cp TOP.nf ${CONDA_PREFIX}/bin/TOP.nf
cp TOP.sh ${CONDA_PREFIX}/bin/TOP
cp nextflow.config ${CONDA_PREFIX}/bin/nexflow.config

conda deactivate

echo "Downloading Kraken database"
echo ""
#wget -O krakenDB.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz
wget -O krakenDB.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20221209.tar.gz

echo "Preparing Kraken database"
echo ""
mkdir krakenDB
tar -xvzf krakenDB.tar.gz -C krakenDB

source activate top_nf
conda env config vars set KRAKENDB=$(pwd)/krakenDB
if [[ ${1} == "--cores" ]]
then
conda env config vars set TOPCORES=${2}

else
conda env config vars set TOPCORES=10

fi

conda deactivate


echo "Downloading docker images"
echo ""
docker pull ghcr.io/garcia-nacho/top_spades
docker pull ghcr.io/garcia-nacho/top_abricate
docker pull ghcr.io/garcia-nacho/top_emmtyper
docker pull ghcr.io/garcia-nacho/top_hicap
docker pull ghcr.io/garcia-nacho/top_seroba
docker pull ghcr.io/garcia-nacho/top_virfinder
docker pull ghcr.io/garcia-nacho/top_prokka

