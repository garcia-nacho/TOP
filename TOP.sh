#!/bin/bash

source activate top_nf

wget -O ${CONDA_PREFIX}/bin/TOP.nf https://github.com/garcia-nacho/TOP/raw/master/TOP.nf
wget -O ${CONDA_PREFIX}/bin/nextflow.config https://github.com/garcia-nacho/TOP/raw/master/nextflow.config

nextflow ${CONDA_PREFIX}/bin/TOP.nf --readsfolder "${1}" --krakenDB ${KRAKENDB}

rm -rf ./work

conda deactivate
