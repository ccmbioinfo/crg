#!/bin/bash
#PBS -l walltime=5:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g,mem=21g

#arguments passed: family_id and sample 

echo "sample = $sample";
if [ -d $sample ]
then
	cd $sample/${family_id}/input;
	ln -s ../../../../remove_decoys/${sample}.no_decoy_reads.bam $sample.bam; 
	cd ../../..;
else
	echo "Decoy file not found"
fi


if [ -d $sample ]
then
	cd $sample;
	echo "current directory: `pwd`";
	~/crg/crg.prepare_bcbio_run.sh ${family_id} sv;
	cd ..;
else
	echo "Participant Folder Not Found."
fi

