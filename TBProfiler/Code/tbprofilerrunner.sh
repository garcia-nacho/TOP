#!/bin/bash
cd /Data
sample=${samplename}
r1=$(ls *_R1*.fastq.gz)
r2=$(ls *_R2*.fastq.gz)
source activate tb-profiler
tb-profiler profile -1 ${r1} -2 ${r2} -p ${sample} --csv
mv results/${sample}.results.csv ${sample}.tb_profiler.csv
mv results/${sample}.results.json ${sample}.tb_profiler.json
rm -rf bam
rm -rf results  
conda deactivate
Rscript /home/docker/Code/TBprofilerparser.R
mv vcf/* .
rm -rf vcf
