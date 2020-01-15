#!/bin/bash

# crg.sv.sh - generates a structural variant report from multiple vcfs
# Dennis Kao, 2019-05-09

# Input: 
#   vcfs = comma seperated paths to vcfs
#   panel = .bed containing regions of interest
#   cohort = a name for the cohort

# Output:
#   .tsv a single report containing annotated information about structural variants across all samples

# Implementation:
#   1. crg.sv.prioritize.sh - filter out SV's using a panel and annotate the VCF using various tools
#   2. crg.intersect_sv_vcfs.sh - assess VCF's for "equivalent" SV's. Then generate an annotated report.

vcfs=$1
panel=$2
cohort=$3

PRIORITIZE_JOB_IDS=()

for vcf in `echo $vcfs | tr ',' ' '`; do 
    PRIORITIZE_JOB_IDS+="afterok:$(qsub ~/crg/crg.sv.prioritize.sh -v vcf=$vcf,panel=$panel),";
    #PRIORITIZE_JOB_IDS+="afterok:$(qsub ./test.sh),";
done

# convert job ids to string
DEPENDENCIES=`echo $PRIORITIZE_JOB_IDS`
DEPENDENCIES="${DEPENDENCIES%?}"

#ANN_VCFS=`ls *svscore*.vcf | tr '\n' ' '`

INTERSECT_JOB_ID=`qsub -F "$cohort" ~/crg/crg.intersect_sv_vcfs.sh -W depend=${DEPENDENCIES}`
