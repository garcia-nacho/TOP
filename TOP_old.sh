#!/bin/bash

source activate top_nf


if [[ ${1} == "--update" ]]
then
    echo "Updating The One Pipeline"
    wget -O ${CONDA_PREFIX}/bin/TOP.nf https://github.com/garcia-nacho/TOP/raw/master/TOP.nf
    wget -O ${CONDA_PREFIX}/bin/nextflow.config https://github.com/garcia-nacho/TOP/raw/master/nextflow.config
    wget -O ${CONDA_PREFIX}/top_template.html https://github.com/garcia-nacho/TOP/raw/master/top_template.html
    wget -O ${CONDA_PREFIX}/bin/TOP.sh https://github.com/garcia-nacho/TOP/raw/master/TOP.sh
    
    docker pull ghcr.io/garcia-nacho/top_spades
    docker pull ghcr.io/garcia-nacho/top_abricate
    docker pull ghcr.io/garcia-nacho/top_emmtyper
    docker pull ghcr.io/garcia-nacho/top_hicap
    docker pull ghcr.io/garcia-nacho/top_seroba
    docker pull ghcr.io/garcia-nacho/top_virfinder
    #docker pull ghcr.io/garcia-nacho/top_prokka
    docker pull ghcr.io/garcia-nacho/top_ngstar
    docker pull ghcr.io/garcia-nacho/top_tbpipeline
    docker pull ghcr.io/garcia-nacho/top_seqsero
    docker pull ghcr.io/garcia-nacho/top_ngmaster
    docker pull ghcr.io/garcia-nacho/top_ecoli
    docker pull ghcr.io/garcia-nacho/top_meningotype
    docker pull ghcr.io/garcia-nacho/top_tartrate

else

    echo "Running The One Pipeline"
    nextflow ${CONDA_PREFIX}/bin/TOP.nf --readsfolder "${1}" --krakenDB "${KRAKENDB}" --TBDB "${TBDB}" --tempfolder "${TEMPDB}" --spadescores ${SPADESCORES} --threads ${TOPCORES} -resume -with-timeline -with-report
    #Delete working directory if there is no error

    if test -f "${1}/TOPresults/Summaries_"*".xlsx"
    then
        echo "Cleaning up..."
        rm -rf work
    fi

    if test -f "${1}/TOPresults/TB_Pipeline/Non_MTBC_samples_in_the_run"
    then
        rm -rf ${1}/TOPresults/TB_Pipeline/
    fi

fi

conda deactivate
