#!/bin/bash

#If ${1} doesnt exist go to pwd
echo "Downloading Kraken database"
echo ""
wget -O ${1}/krakenDB https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz

echo "Unpacking Kraken database"
echo ""
tar -xvzf ${1}/krakenDB https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220926.tar.gz

echo "Installing docker images"
echo ""
docker build -t garcianacho/seroba ./Seroba/
docker build -t garcianacho/prokka ./Prokka/
docker build -t garcianacho/top ./Spades/


