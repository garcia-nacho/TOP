#!/bin/bash
Tab=$(printf '\t')
source activate hicap

for files in *.fasta
do

hicap -q ${files} -o Results --threads 8 --log_fp ${files%.fasta}.log
if ! test -f "Results/${files%.fa*}.tsv"; then
    echo "#isolate${Tab}predicted_serotype${Tab}attributes${Tab}genes_identified${Tab}locus_location${Tab}region_I_genes${Tab}region_II_genes${Tab}region_III_genes${Tab}IS1016_hits" > Results/${files%.fa*}.tsv
    echo "${files%.fa*}${Tab}NTHi${Tab}NA${Tab}NA${Tab}NA${Tab}NA${Tab}NA${Tab}NA" >> Results/${files%.fa*}.tsv
fi

done

cp Results/* ./

conda deactivate
