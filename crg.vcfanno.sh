#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $vcf ]
then
    vcf=$1
fi

bname=`basename $vcf .vcf.gz`

prefix=$HOME/crg

vcfanno -p 5 -lua $prefix/crg.vcfanno.lua \
	     -base-path /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/variation \
	     $prefix/crg.vcfanno.conf \
	     $vcf | sed -e 's/Number=A/Number=1/g' | bgzip -c > $bname.annotated.vcf.gz

tabix $bname.annotated.vcf.gz
