#!/bin/bash

#Usage: bash downloadHG19.sh PATH_TO_DIR_TO_PLACE_REFERENCE_GENOME
set -exu
HG38_DIR=$1
mkdir -p $HG38_DIR
cd $HG38_DIR

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar -zxvf hg38.chromFa.tar.gz
for i in chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrY.fa chrX.fa chrM.fa; do cat chroms/$i >> hg38.fa ; rm chroms/$i ; done
rm -rf chroms
samtools faidx hg38.fa

#in config.yaml update the path to the downloaded HG38.fa
