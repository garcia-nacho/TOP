#!/bin/bash

echo Downloading cgMLST profiles from pasteur.fr
GET https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/4/profiles_csv > profiles_cgmlst.tsv

for files in $(ls *.fa*)
do

echo Analyzing ${files%.fa*}
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/4/sequence" -d @- > ${files%.fa*}_cgmlst.json
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3/sequence" -d @- > ${files%.fa*}_mlst.json

done

echo Putting things together
Rscript /home/docker/code/mlstparser.R
