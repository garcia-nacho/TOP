#!/bin/sh

source activate top_nf

#Chech if ${1} exists otherwise 1=$(pwd)
TOP.nf -readsfolder "${1}"
rm -rf ./work

conda deactivate