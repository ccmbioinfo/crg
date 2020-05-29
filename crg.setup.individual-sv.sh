#!/bin/bash
#PBS -l walltime=5:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g,mem=21g

for f in *; do cd $f/${family_id}/input;
	if [ -d $f ]
	then 
		ln -s ../../../../remove_decoys/$f.no_decoy_reads.bam $f.bam; 
		cd ../../..;
	else
		echo "decoy file not found"
	fi

done

for f in *; 
do 
	if [ -d $f ]
	then
		cd $f; 
		~/crg/crg.prepare_bcbio_run.sh ${family_id} sv;
		cd ..;
	else
		echo "Participant Folder Not Found."
	fi
done