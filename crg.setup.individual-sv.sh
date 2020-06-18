#!/bin/bash
#PBS -l walltime=5:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g,mem=21g

for f in ${family_id}*; 
do
	echo $f
	if [ -d $f ]
	then 
		cd $f/${family_id}/input;
		ln -s ../../../../remove_decoys/${f}-ready.no_decoy_reads.bam $f.bam; 
		cd ../../..;
	else
		echo "Decoy file not found"
	fi

done

for f in ${family_id}*; 
do 
	if [ -d $f ]
	then
		cd $f;
		pwd
		~/crg/crg.prepare_bcbio_run.sh ${family_id} sv;
		cd ..;
	else
		echo "Participant Folder Not Found."
	fi
done