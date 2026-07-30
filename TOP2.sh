#!/bin/bash

source activate top2_nf
SHORT=u,r:,f:,h,t:,c:,d,s,x,l,n
LONG=update,reads:,fastas:,help,tbdb:,cores:,dev,uninstall,clean,soft,ont
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
    echo "-s or --soft to run a soft trimming of the adapters"
    echo "-o or --ont to run the Nanopore version of the pipeline"
    echo ""
    conda deactivate
    exit 1
fi

SOFT="No"
ONT="No"
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
    -t | --tbdb )
      TBDB="$2"
      shift 2
      ;;
    -r | --reads )
      READS="$2"
      shift 2
      ;;
    -s | --soft )
      SOFT="Yes"
      shift
      ;;
    -o | --ont )
      ONT="Yes"
      shift
      ;;
    -x | --uninstall )
      echo "Removing TOP2"
      rm ${TOPSHPATH}/TOP2.sh
      conda deactivate top2_nf
      conda remove -n top2_nf --all
      echo "Note that you will have to delete the docker containers manually"
      exit 0
      ;;
    -f | --fastas )
      echo "Running from fastas is not implemented yet"
      exit 0
      ;;
    -l | --clean )
      
      echo "Cleaning TOP2 temporary files"

      rm -rf work
      rm .nextflow.lo*
      rm -rf .nexflo*
      rm -rf $(pwd)/top_progress
      rm report-*.html
      rm timeline-*.html

      exit 0

      ;;
    -u | --update )
      echo "Updating The One Pipeline"
      wget -O ${CONDA_PREFIX}/bin/TOP.nf https://github.com/garcia-nacho/TOP/raw/dev/TOP.nf
      wget -O ${CONDA_PREFIX}/bin/TOP_ont.nf https://github.com/garcia-nacho/TOP/raw/dev/TOP_ont.nf
      wget -O ${CONDA_PREFIX}/bin/nextflow.config https://github.com/garcia-nacho/TOP/raw/dev/nextflow.config
      wget -O ${CONDA_PREFIX}/top_template.html https://github.com/garcia-nacho/TOP/raw/dev/top_template.html
      wget -O ${CONDA_PREFIX}/bin/TOP2.sh https://github.com/garcia-nacho/TOP/raw/dev/TOP2.sh
      
      docker pull ghcr.io/garcia-nacho/top_spades:v.1.2
      docker pull ghcr.io/garcia-nacho/top_abricate:v1.1
      docker pull ghcr.io/garcia-nacho/top_emmtyper:v1.1
      docker pull ghcr.io/garcia-nacho/top_hicap:v1.1
      #docker pull ghcr.io/garcia-nacho/top_seroba:v1.1
      docker pull ghcr.io/garcia-nacho/top_virfinder:v1.1
      #docker pull ghcr.io/garcia-nacho/top_prokka
      docker pull ghcr.io/garcia-nacho/top_ont
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
      docker pull ghcr.io/garcia-nacho/top_pbp3
      docker pull ghcr.io/garcia-nacho/top_ont
      docker pull ghcr.io/garcia-nacho/top_seroba2
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
      echo "use -n or --ont to run the Nanopore version of the pipeline"
      echo "use -r or --reads to run from a set of paired fastq files e.g. TOP2.sh -f /path/to/fastq"
      echo "use -c or --cores to set the number of cores"
      echo "use -t or --tbdb to set the path to the mtb database"
      echo "use -f or --fastas to run from a set of fastas e.g. TOP.sh -f /path/to/fastas"
      echo "use -d or --dev to run in development mode.(i.e all files are saved after a successful run)" 
      echo "use -h or --help to show this help" 
      echo "use -u or --update to update the pipeline"
      echo "use -p or --private to prevent updating the database on test runs or private samples"
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

mkdir top_progress

if ${ONT}=="No"
then

  echo "Running The One Pipeline V2: Illumina"
  echo ""
  echo nextflow ${CONDA_PREFIX}/bin/TOP.nf --readsfolder "${READS}" --krakenDB "${KRAKENDB}" --TBDB "${TBDB}" --tempfolder "${TEMPDB}" --spadescores ${SPADESCORES} --threads ${TOPCORES} --softtrimming ${SOFT} -resume -with-timeline -with-report
  echo ""
  nextflow ${CONDA_PREFIX}/bin/TOP.nf --readsfolder "${READS}" --krakenDB "${KRAKENDB}" --TBDB "${TBDB}" --tempfolder "${TEMPDB}" --spadescores ${SPADESCORES} --threads ${TOPCORES} --progress_dir $(pwd)/top_progress --softtrimming ${SOFT} -resume -with-timeline -with-report

else

  echo "Running The One Pipeline V2: Nanopore"
  echo ""
  echo nextflow ${CONDA_PREFIX}/bin/TOP_ont.nf --readsfolder "${READS}" --krakenDB "${KRAKENDB}" --TBDB "${TBDB}" --tempfolder "${TEMPDB}" --spadescores ${SPADESCORES} --threads ${TOPCORES} --softtrimming ${SOFT} -resume -with-timeline -with-report
  echo ""

  nextflow ${CONDA_PREFIX}/bin/TOP_ont.nf --readsfolder "${READS}" --krakenDB "${KRAKENDB}" --TBDB "${TBDB}" --tempfolder "${TEMPDB}" --spadescores ${SPADESCORES} --threads ${TOPCORES} --progress_dir $(pwd)/top_progress --softtrimming ${SOFT} -resume -with-timeline -with-report

fi


if test -f "${READS}/TOPresults/Summaries_"*".xlsx"
then
    echo "Cleaning up..."
    rm -rf work
    rm -rf top_progress
fi

if test -f "${READS}/TOPresults/TB_Pipeline/Non_MTBC_samples_in_the_run"
then
    rm -rf ${READS}/TOPresults/TB_Pipeline/
fi

conda deactivate