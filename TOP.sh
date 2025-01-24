#!/bin/bash

source activate top_nf
SHORT=u,r:,f:,h,t:,c:,d,s:,x
LONG=update,reads:,fastas:,help,tbdb:,cores:,dev,sp:,uninstall
OPTS=$(getopt --options $SHORT --longoptions $LONG -- "$@")

if [ $? != 0 ]; then
    echo ""
    echo "Error with the option/s provided" >&2
    echo "use -r or --reads to run from a set of paired fastq files e.g. TOP.sh -f /path/to/fastq"
    echo "use -c or --cores to set the number of cores"
    echo "use -t or --tbdb to set the path to the mtb database"
    echo "use -f or --fastas to run from a set of fastas e.g. TOP.sh -f /path/to/fastas"
    echo "use -d or --dev to run in development mode.(i.e all files are saved after a successful run)" 
    echo "use -h or --help to show this help" 
    echo "use -u or --update to update the pipeline"
    echo "use -x or --uninstall to uninstall the pipeline" 
    echo ""
    conda deactivate top_nf
    exit 1
fi


READS=$(pwd)
devmode=0
eval set -- "$OPTS"

while :
do
  case "$1" in
    -c | --cores )
      TOPCORES="$2"
      SPADESCORES=$((${TOPCORES}-2))
      shift 2
      ;;
    -s | --sp )
      echo "Specie specific tools are not implemented yet"
      exit 2
      ;;
    -t | --tbdb )
      TBDB="$2"
      shift 2
      ;;
    -r | --reads )
      READS="$2"
      shift 2
      ;;
    -x | --uninstall )
      echo "Removing TOP"
      rm ${TOPSHPATH}/TOP.sh
      conda deactivate top_nf
      conda remove -n top_nf --all
      echo "Note that you will have to delete the docker containers manually"
      exit 0
      ;;
    -f | --fastas )
      echo "Running from fastas is not implemented yet"
      exit 0
      ;;
    -u | --update )
      echo "Updating The One Pipeline"
      wget -O ${CONDA_PREFIX}/bin/TOP.nf https://github.com/garcia-nacho/TOP/raw/master/TOP.nf
      wget -O ${CONDA_PREFIX}/bin/nextflow.config https://github.com/garcia-nacho/TOP/raw/master/nextflow.config
      wget -O ${CONDA_PREFIX}/top_template.html https://github.com/garcia-nacho/TOP/raw/master/top_template.html
      wget -O ${CONDA_PREFIX}/bin/TOP.sh https://github.com/garcia-nacho/TOP/raw/master/TOP.sh
      
      docker pull ghcr.io/garcia-nacho/top_spades:v.1.1
      docker pull ghcr.io/garcia-nacho/top_abricate:v1.1
      docker pull ghcr.io/garcia-nacho/top_emmtyper:v1.1
      docker pull ghcr.io/garcia-nacho/top_hicap:v1.1
      docker pull ghcr.io/garcia-nacho/top_seroba:v1.1
      docker pull ghcr.io/garcia-nacho/top_virfinder:v1.1
      #docker pull ghcr.io/garcia-nacho/top_prokka
      docker pull ghcr.io/garcia-nacho/top_ngstar:v1.1
      docker pull ghcr.io/garcia-nacho/top_tbpipeline:v1.1
      docker pull ghcr.io/garcia-nacho/top_seqsero:v1.1
      docker pull ghcr.io/garcia-nacho/top_ngmaster:v1.1
      docker pull ghcr.io/garcia-nacho/top_ecoli:v1.1
      docker pull ghcr.io/garcia-nacho/top_meningotype:v1.1
      docker pull ghcr.io/garcia-nacho/top_tartrate:v1.1
      docker pull ghcr.io/garcia-nacho/top_amrfinderplus
      docker pull ghcr.io/garcia-nacho/top_tbprofiler
      docker pull ghcr.io/garcia-nacho/top_diphtoscan
      docker pull ghcr.io/garcia-nacho/top_bpprofiler


      exit 0
      ;;
    -d | --dev )
      devmode=1
      shift
      ;;
    -h | --help)
      echo ""
      echo "This is TOP (The One Pipeline) "
      echo ""
      echo "use -r or --reads to run from a set of paired fastq files e.g. TOP.sh -f /path/to/fastq"
      echo "use -c or --cores to set the number of cores"
      echo "use -t or --tbdb to set the path to the mtb database"
      echo "use -f or --fastas to run from a set of fastas e.g. TOP.sh -f /path/to/fastas"
      echo "use -d or --dev to run in development mode.(i.e all files are saved after a successful run)" 
      echo "use -h or --help to show this help" 
      echo "use -u or --update to update the pipeline"
      echo "use -x or --uninstall to uninstall the pipeline" 
      echo ""
      exit 0
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      exit 2
      ;;
  esac
done


echo "Running The One Pipeline V1.1:"
echo ""
echo nextflow ${CONDA_PREFIX}/bin/TOP.nf --readsfolder "${READS}" --krakenDB "${KRAKENDB}" --TBDB "${TBDB}" --tempfolder "${TEMPDB}" --spadescores ${SPADESCORES} --threads ${TOPCORES} -resume -with-timeline -with-report
#Delete working directory if there is no error
echo ""
nextflow ${CONDA_PREFIX}/bin/TOP.nf --readsfolder "${READS}" --krakenDB "${KRAKENDB}" --TBDB "${TBDB}" --tempfolder "${TEMPDB}" --spadescores ${SPADESCORES} --threads ${TOPCORES} -resume -with-timeline -with-report

if test -f "${READS}/TOPresults/Summaries_"*".xlsx"
then
    if devmode == 0
    then
    echo "Cleaning up..."
    rm -rf work
    fi
fi

if test -f "${READS}/TOPresults/TB_Pipeline/Non_MTBC_samples_in_the_run"
then
    rm -rf ${READS}/TOPresults/TB_Pipeline/
fi

conda deactivate

