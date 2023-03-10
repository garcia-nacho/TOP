#!/bin/sh

source activate top_nf

nextflow ${CONDA_PREFIX}/bin/TOP.nf --readsfolder "${1}" --krakenDB ${CONDA_PREFIX}/krakenDB
rm -rf ./work

conda deactivate