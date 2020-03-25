#!/bin/bash

#PBS -l walltime=0:20:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

IN=$1
OUT="${IN}.snpeff.vcf"

/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/snpEff -Xms750m -Xmx20g  -i vcf -o vcf -dataDir /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/snpeff \
GRCh37.75 \
${IN} > ${OUT}
