#!/bin/bash

#for each of the ERDS calls, create a bed file --> temp.bed
#look for bed intersection with "-a" --> proband "-b" set of family bed

#for each entry in a, write out entry in A + overlapping entries in b/c by name
#this file will be used to compare each A entry and add B/C entries to it

sample=$1
proband=$2
f1=$3
f2=$4

#run: sh merge.ERDS_Cnvnator.with.bedtools.sh 1159R 1159R_CH0531 1159R_CH0880 1159R_CH0879 

ls /Users/arun/Documents/Projects/C4R/WGS.SV/$sample/*erds*tsv > $sample.erds_cnvnator.list

cat $sample.erds_cnvnator.list | while read -r file
do
    echo $file
    cat $file | awk {'printf ("%s\t%s\t%s\t%s\n", $2,$3,$4,$1)'} | grep -v "CHROM" > $file.bed
done

cd /Users/arun/Documents/Projects/C4R/WGS.SV/$sample

bedtools intersect -a $proband.erds\+_cnvAnn_rev1.0.0.tagged.tsv.bed -b $f1.erds\+_cnvAnn_rev1.0.0.tagged.tsv.bed $f2.erds\+_cnvAnn_rev1.0.0.tagged.tsv.bed -names $f1 $f2  -wa -wb -f 0.5 > $sample.intersect.bed



