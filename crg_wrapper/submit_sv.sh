#!/bin/bash

#PBS -N submit_sv
#PBS -l nodes=1:ppn=1
#PBS -joe 
#PBS -d .
#PBS -l vmem=20g,mem=20g
#PBS -m ae


#submit this from base project directory: /hpf/largeprojects/ccmbio/aarthi/proj_CHEO/CRG/496
#usage: qsub submit_sv.sh -F "496" 

family_id=$1;
project_dir=`pwd`; #base project path
logfile="${project_dir}/bcbio-sv/${family_id}_sv_jobids.log";

if [ -f $logfile ]; then	
	rm $logfile
fi;
touch "$logfile"


#step1: remove decoys and call svs 
#run this from base project directory
~/crg/crg.call-svs.sh ${family_id}

#read jobids written by above script in logfile
if [ -f "$logfile" ]; then
	declare -A jobarr;
	while IFS== read -r key value; do
		jobarr[$key]=$value;
	done < "$logfile"
fi
echo "sv bcbio jobids: ${jobarr["sv_bcbio"]}"

#step2: merge, annotate and intersect: this requires project spec HPO file
#if not found, print the message and still run 
#run this from bcbio-sv
cd bcbio-sv
echo "Entering dir=`pwd`";
if [ ! -f  "~/gene_data/HPO/${family_id}_HPO.txt" ]; then 
	echo "${family_id}_HPO.txt file not found in ~/gene_data/HPO folder, crg.merge.annotate.sv.sh will still be run";
	#~/crg/crg.merge.annotate.sv.sh script won't be run";
fi;
qsub ~/crg/crg.merge.annotate.sv.sh  -F "${family_id}" -W depend=afterok:${jobarr[sv_bcbio]}

#step3: prioritize 
#this requires panel, panel100kflank from small-variants calling 
#this must also be run after jobids from "~/crg/crg.call-svs.sh 1582"  finishes
#and creates required directory for looping
#best to put this in seperate script with depend directive

#for i in ${project_dir}"/"${family_id}"/bcbio-sv/"${family_id}_*; do
#	cd $i;
#	echo "entering: `pwd`";
#	#qsub ~/crg/crg.sv.prioritize.sh -v case=${family_id},panel=../${family_id}.bed 
#	cd $cwd
#	echo "exiting to: `pwd`"
#done;

