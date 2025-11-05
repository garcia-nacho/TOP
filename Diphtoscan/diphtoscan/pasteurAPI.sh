echo Downloading profiles from pasteur.fr

GET https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef/schemes/3/profiles_csv > profiles_mlst.tsv


for files in $(ls *.fasta)
do

echo Analyzing ${files%.fasta}
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef/schemes/3/sequence" -d @- > ${files%.fasta}_mlst.json
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef/schemes/4/sequence" -d @- > ${files%.fasta}_toxinst.json
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef/schemes/6/sequence" -d @- > ${files%.fasta}_spuA_island.json
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef/schemes/2/sequence" -d @- > ${files%.fasta}_pbp.json

done

echo Putting things together

