#!/bin/bash
Tab=$(printf '\t')
source activate hicap
mkdir Results
for files in *.fasta
do

hicap -q ${files} -o Results --threads 8 --log_fp ${files%.fasta}.log
if ! test -f "Results/${files%.fa*}.tsv"; then
    echo "#isolate${Tab}predicted_serotype${Tab}attributes${Tab}genes_identified${Tab}locus_location${Tab}region_I_genes${Tab}region_II_genes${Tab}region_III_genes${Tab}IS1016_hits" > Results/${files%.fa*}.tsv
    echo "${files%.fa*}${Tab}NTHi${Tab}NA${Tab}NA${Tab}NA${Tab}NA${Tab}NA${Tab}NA" >> Results/${files%.fa*}.tsv
fi

done

if test -f "Results/${files%.fa*}.gbk"; then mv  Results/${files%.fa*}.gbk ./ Results/${files%.fa*}_HiCap.gbk; fi
if test -f "Results/${files%.fa*}.tsv"; then mv  Results/${files%.fa*}.tsv ./ Results/${files%.fa*}_HiCap.tsv; fi
if test -f "Results/${files%.fa*}.log"; then mv  Results/${files%.fa*}.tsv ./ Results/${files%.fa*}_HiCap.log; fi
if test -f "Results/${files%.fa*}.svg"; then mv  Results/${files%.fa*}.svg ./ Results/${files%.fa*}_HiCap.svg; fi

rm -rf Results
conda deactivate
