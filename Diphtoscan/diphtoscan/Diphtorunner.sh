#!/bin/bash
mkdir diphtoresults
for files in $(ls *.fa*)
do
echo ${files%.fa*} 

python /home/docker/diphtoscan/__main__.py -a ${files} --mlst --tox -res_vir -plus -o diphtoresults/${files%.fa*} 

done

cd diphtoresults
Rscript /home/docker/diphtoscan/diphtoparser.R
