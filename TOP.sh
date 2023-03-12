#!/bin/bash

source activate top_nf

if [[ ${1} == "--update" ]]
then
echo "Updating The One Pipeline"
wget -O ${CONDA_PREFIX}/bin/TOP.nf https://github.com/garcia-nacho/TOP/raw/master/TOP.nf
wget -O ${CONDA_PREFIX}/bin/nextflow.config https://github.com/garcia-nacho/TOP/raw/master/nextflow.config
docker pull ghcr.io/garcia-nacho/top_spades
docker pull ghcr.io/garcia-nacho/top_abricate
docker pull ghcr.io/garcia-nacho/top_emmtyper
docker pull ghcr.io/garcia-nacho/top_hicap
docker pull ghcr.io/garcia-nacho/top_seroba
docker pull ghcr.io/garcia-nacho/top_virfinder
docker pull ghcr.io/garcia-nacho/top_prokka
else
echo "Running The One Pipeline"
wget -O ${CONDA_PREFIX}/bin/nextflow.config https://github.com/garcia-nacho/TOP/raw/master/nextflow.config
nextflow ${CONDA_PREFIX}/bin/TOP.nf --readsfolder "${1}" --krakenDB ${KRAKENDB}
rm -rf ./work
fi

conda deactivate
