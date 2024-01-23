#!/bin/bash
mkdir /results

for dir in $(ls -d */)
do
cd ${dir}

R1=$(ls *_R1*.fastq.gz)
R2=$(ls *_R2*.fastq.gz)
SeqSero2_package.py -p 10 -t 2 -i ${R1} ${R2}
mv $(ls -d */) /results/${dir%/}_SeqSeroResults_allele

SeqSero2_package.py -m k -t 2 -i ${R1} ${R2}
mv $(ls -d */) /results/${dir%/}_SeqSeroResults_kmer

cd ..
done

cp -r /results ./SeqSero2.results

cd SeqSero2.results
Rscript /Code/postanalysis.R

