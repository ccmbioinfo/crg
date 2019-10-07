#!/bin/bash
#PBS -l walltime=5:00:00,nodes=1:ppn=2
#PBS -joe .
#PBS -d .
#PBS -l vmem=8g,mem=8g

#Given a set of vcfs, merge them with SURVIVOR
#
#iterate through each length bin and apply
#these are in this file survivor_merge_bins.txt
#min max distance
#0 100 50
#100 300 150
#300 500 250
#500 1000 500
#1000 2500 1250
#2500 5000 2500
#5000 10000 5000
#10000 -1 6000
#
#
#

#combine vcf files to create a vcf.list

SAMPLE=$1
MERGED_SAMPLE="merged.${SAMPLE}.vcf"
TMP_VCF="temp.vcf"
TMP_VCF_LIST="tmp.vcf.list"
BINS="/home/dennis.kao/crg/vcf_merging/survivor_merge_bins.txt"

#~/Apps/SURVIVOR/Debug/SURVIVOR filter
#VCF file to filter
#BED file with regions to ignore (NA to disable)
#Min SV size (-1 to disable)
#Max SV size (-1 to disable)
#Min allele frequency (0-1)
#Min number of reads support: RE flag (-1 to disable)
#Output vcf file

#~/Apps/SURVIVOR/Debug/SURVIVOR merge
#File with VCF names and paths
#max distance between breakpoints 
#Minimum number of supporting caller
#Take the type into account (1==yes, else no)
#Take the strands of SVs into account (1==yes, else no)
#Estimate distance based on the size of SV (1==yes, else no).
#Minimum size of SVs to be taken into account.
#Output VCF filename


#a --> min
#b --> max
#c --> distance

if [ ! -e ${MERGED_SAMPLE} ]
then
    rm ${MERGED_SAMPLE}
fi


cat ${BINS} | while read -r a b c
do
    for vcf in *-metasv.vcf; do
        echo ${vcf}
        SURVIVOR filter ${vcf} NA $a $b 0 -1 ${vcf}.tmpvcf
    done

    ls ./*.tmpvcf > ${TMP_VCF_LIST}
    SURVIVOR merge ${TMP_VCF_LIST} ${c} 1 1 1 1 -1 ${TMP_VCF}

    if [ ! -e ${MERGED_SAMPLE} ]
    then
        grep "#" ${TMP_VCF} > ${MERGED_SAMPLE}
    fi
    grep -v "#" ${TMP_VCF} >> ${MERGED_SAMPLE}
done

rm ${TMP_VCF}
rm *.tmpvcf