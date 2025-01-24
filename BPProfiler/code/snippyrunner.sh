#!/bin/bash
mkdir /Data/SnippyResults

ref=/Data/PertussisRoI.fa

for dir in $(ls -d */)
do

     cd ${dir}
     echo ${dir}
        
     R1=$(ls *R1*)
     R2=$(ls *R2*)
        
     snippy --cpus 6 --outdir /Data/SnippyResults/${dir%/} --ref ${ref} --R1 ${R1} --R2 ${R2}
     cd ..

done
cd /Data/SnippyResults/
#java -Xmx8g -jar /home/docker/snpEff/snpEff.jar -v Bordetella_pertussis_tohama_i  SnippyResults/2305704-BO-24GB002696/snps.filt.vcf > test_ann.vcf

