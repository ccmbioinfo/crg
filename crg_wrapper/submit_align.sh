#!/bin/bash


#PBS -N submit_align
#PBS -joe
#PBS -d .
#PBS -l nodes=1:ppn=1
#PBS -l vmem=20g,mem=20g
#PBS -m ae

#currently not using this script in crg_wrapper.py
#submit this script from project/bcbio-align directory
#usage: qsub submit_align.sh -F "family analysistype"


read -r family analysis <<< $(echo $@ | awk '{print $1, $2;}')
family_id=$family;

if [ -z $analysis ]; then 
    analysis="align_decoy";
fi;

wd=`pwd`;
logfile="${wd}/${family_id}_jobids.log"; #all align job ids are recorded here
touch ${logfile};

#step1: setup bcbio directories for align
echo "Setting directories and config required for align with decoy: ${family_id}";
~/crg/crg.prepare_bcbio_run.sh ${family_id} $analysis

#step2: submit align with decoy job
bcbio_job=$(qsub ~/cre/bcbio.pbs -v project=${family_id});
echo "submitted align with decoy bcbio jobid for ${family_id}: ${bcbio_job}";

#record submitted jobids here
echo "align_bcbio=${bcbio_job}" >> ${logfile};

