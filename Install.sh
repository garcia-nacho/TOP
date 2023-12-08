#!/bin/bash

SHORT=t:,c:,h
LONG=tbdb:,cores:,help
OPTS=$(getopt --options $SHORT --longoptions $LONG -- "$@")

eval set -- "$OPTS"

while :
do
  case "$1" in
    -c | --cores )
      cores="$2"
      shift 2
      ;;
    -t | --tbdb )
      tbdb="$2"
      shift 2
      ;;
    -k | --kraken )
      kraken="$2"
      shift 2
      ;;
    -h | --help)
      echo "This is a TOP installation script"
      echo "use -c or --cores to set the number of cores"
      echo "use -t or --tbdb to set the path to the tb_database"
      echo "use -k or --kraken to set the path to the tb_database"
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

#Cores
if [ -z ${cores+x} ]
then
${cores}=10
fi
echo "Number of cores assigned to TOP: " ${cores}
echo ""

#TBDB
if [ -z ${tbdb+x} ]
then
mkdir $(pwd)/TBDB
${tbdb}=$(pwd)/TBDB
echo "Creating TBDB in "$(pwd)/TBDB
echo ""
fi

conda create -n top_nf -y -c bioconda nextflow
source activate top_nf
cp TOP.nf ${CONDA_PREFIX}/bin/TOP.nf
cp TOP.sh ${CONDA_PREFIX}/bin/TOP
cp nextflow.config ${CONDA_PREFIX}/bin/nexflow.config
mkdir ${CONDA_PREFIX}/top_temp
conda deactivate

if [ -z ${kraken+x} ]
then
echo "No Kraken database set... Downloading it into"$(pwd)/krakenDB
echo ""
wget -O krakenDB.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20221209.tar.gz
mkdir krakenDB
echo "Unpacking Kraken database"
echo ""
tar -xvzf krakenDB.tar.gz -C krakenDB
kraken=$(pwd)/krakenDB
fi


source activate top_nf
conda env config vars set KRAKENDB=${kraken}
conda env config vars set TBDB=${tbdb}
conda env config vars set TempDB=${CONDA_PREFIX}/top_temp
conda env config vars set TOPCORES=${2}

conda deactivate

echo "Downloading docker images"
echo ""
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
docker pull push ghcr.io/garcia-nacho/top_tbpipeline
