#!/bin/bash

source activate top_nf
SHORT=u,r:,f:,h,t:,c:,d
LONG=update,reads:,fastas:,help,tbdb:,cores:,dev
OPTS=$(getopt --options $SHORT --longoptions $LONG -- "$@")

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
    -f | --fastas )
      echo "Running from fastas is not implemented yet"
      exit 2
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
      echo ""
      exit 2
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      ;;
  esac
done


echo "Running The One Pipeline:"
echo ""
echo nextflow ${CONDA_PREFIX}/bin/TOP.nf --readsfolder "${READS}" --krakenDB "${KRAKENDB}" --TBDB "${TBDB}" --tempfolder "${TEMPDB}" --spadescores ${SPADESCORES} --threads ${TOPCORES} -resume -with-timeline -with-report
#Delete working directory if there is no error
echo ""
nextflow ${CONDA_PREFIX}/bin/TOP.nf --readsfolder "${READS}" --krakenDB "${KRAKENDB}" --TBDB "${TBDB}" --tempfolder "${TEMPDB}" --spadescores ${SPADESCORES} --threads ${TOPCORES} -resume -with-timeline -with-report

if test -f "${READS}/TOPresults/Summaries_"*".xlsx"
then
    echo "Cleaning up..."
    rm -rf work
fi

if test -f "${READS}/TOPresults/TB_Pipeline/Non_MTBC_samples_in_the_run"
then
    rm -rf ${READS}/TOPresults/TB_Pipeline/
fi

conda deactivate

