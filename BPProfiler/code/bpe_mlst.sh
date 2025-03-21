#!/bin/bash

echo Downloading cgMLST profiles from pasteur.fr

GET https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3/profiles_csv > profiles_mlst.tsv
GET https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/4/profiles_csv > profiles_cgmlst.tsv
#GET https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/7/profiles_csv > profiles_BPagST.tsv
#GET https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/10/profiles_csv > profiles_tss3.tsv

for files in $(ls *.fa*)
do

echo Analyzing ${files%.fa*}
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/4/sequence" -d @- > ${files%.fa*}_cgmlst.json
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/3/sequence" -d @- > ${files%.fa*}_mlst.json
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/7/sequence" -d @- > ${files%.fa*}_bpagst.json
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/10/sequence" -d @- > ${files%.fa*}_tss3st.json

(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/8/sequence" -d @- > ${files%.fa*}_phasest.json
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/9/sequence" -d @- > ${files%.fa*}_toxinst.json
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_bordetella_seqdef/schemes/6/sequence" -d @- > ${files%.fa*}_amrst.json


done

echo Putting things together
Rscript /home/docker/code/mlstparser.R
