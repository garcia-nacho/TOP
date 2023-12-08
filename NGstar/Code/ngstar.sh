#!/bin/bash

cd /home/docker/pyngSTar/ 
python3 pyngSTar.py -f -a -i /Data/*.fa* -p pyngSTarDB_050923/ -o ngstar_results.tsv 
mkdir /Data/Seqs23S
cp *.23S.* /Data/Seqs23S
cp ngstar_results.tsv /Data/ngstar_results.tsv
