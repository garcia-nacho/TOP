#!/bin/bash

SHORT=k:,t:,c:,h:,p:
LONG=kraken:,tbdb:,cores:,help:,path:
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
    -p | --path )
      binpath="$2"
      shift 2
      ;;
    -h | --help)
      echo "This is a TOP installation script"
      echo "use -c or --cores to set the number of cores"
      echo "use -t or --tbdb to set the path to the tb_database"
      echo "use -k or --kraken to set the path to the tb_database"
      echo "use -p or --path to set the path to install the main script" 
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
fi
echo "Using TBDB located in "${tbdb}
echo ""

conda create -n top_nf -y
source activate top_nf
conda install -c bioconda nextflow
cp TOP.nf ${CONDA_PREFIX}/bin/TOP.nf

cp nextflow.config ${CONDA_PREFIX}/bin/
cp top_template.html ${CONDA_PREFIX}/top_template.html
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
echo "Using krakenDB located in "${kraken}
echo ""


source activate top_nf
conda env config vars set KRAKENDB=${kraken}
conda env config vars set TBDB=${tbdb}
conda env config vars set TEMPDB=${CONDA_PREFIX}/top_temp
conda env config vars set TOPCORES=${cores}
conda env config vars set SPADESCORES=$((${cores}-2))
if [ -z ${binpath+x} ]
then
cp TOP.sh ${CONDA_PREFIX}/bin/TOP.sh
ln -s ${CONDA_PREFIX}/bin/TOP.sh ${binpath}/TOP.sh
conda env config vars set TOPSHPATH=${binpath}

else
mv TOP.sh ${CONDA_PREFIX}/bin/TOP.sh
ln -s ${CONDA_PREFIX}/bin/TOP.sh $(pwd)/TOP.sh
conda env config vars set TOPSHPATH=$(pwd)

fi

conda deactivate
