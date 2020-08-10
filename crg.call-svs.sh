#!/bin/bash

# input is a family id
# run from the top level folder for the project (i.e. within /<project>) 

family_id=$1
currdir=`pwd`;
logfile="${currdir}/bcbio-sv/${family_id}_sv_jobids.log";

if [ ! -d "bcbio-align" ]
then
	echo "bcbio-align folder not found"
	exit
fi

if [ -d "remove_decoys" ]
then
	echo "Decoys folder already exists. Will skip removing decoys."
else
	echo "Linking bam files from bcbio-align/${family_id} and Submitting Jobs to Remove Decoy Reads"
	
	mkdir remove_decoys
	cd remove_decoys
	
	decoy_jobs=() # capture job ids
	for f in ../bcbio-align/${family_id}/final/${family_id}_*;
	do
		sample=$(basename ${f})
		bam_file="../bcbio-align/${family_id}/final/${sample}/${sample}-ready.bam"
		if [ -f ${bam_file} ];
		then
			ln -s ${bam_file} ${sample}.bam
			decoy_jobs+=($(qsub ~/cre/cre.bam.remove_decoy_reads.sh -v bam=${sample}.bam)) 
		else
			echo "Aligned bam file not found for sample ${bam_file}"
			exit
		fi
	done
	echo "Decoys being removed with job ids: "$( IFS=$', '; echo "${decoy_jobs[*]}" )
	decoy_string=$( IFS=$':'; echo "${decoy_jobs[*]}" ) # store job ids
	
	cd ..
fi

echo "Preparing and submitting SV calling jobs"
cd bcbio-sv
bcbio_jobs=()
for f in ../bcbio-align/${family_id}/final/${family_id}_*; 
do
	sample=$(basename ${f})
	if [ -d ${sample} ];
	then
		echo "Directory structure already exists in bcbio-sv, please remove and re-run"
		exit
	fi	
	mkdir -p ${sample}/${family_id}/input
	if [ -z ${decoy_string} ];
	then
		#pass $family_id and $sample to crg.setup.individual-sv.sh
		setup_bcbio=$(qsub ~/crg/crg.setup.individual-sv.sh -v family_id="${family_id}",sample="${sample}") 
	else
		setup_bcbio=$(qsub ~/crg/crg.setup.individual-sv.sh -v family_id="${family_id}",sample="${sample}" -W depend=afterany:"${decoy_string}")
	fi
	# submit bcbio job to run after setup complete
	cd ${sample}
	bcbio_jobs+=($(qsub ~/cre/bcbio.pbs -v project="${family_id}" -W depend=afterany:"${setup_bcbio}"))
	cd ..
done
bcbio_string=$( IFS=$':'; echo "${bcbio_jobs[*]}" )
echo "SV calling job ids: "$( IFS=$', '; echo "${bcbio_string}" )
echo "Setup Complete."

#write jobids to log file
echo "decoy=${decoy_string}" >> $logfile
echo "sv_setup=${setup_bcbio}" >> $logfile
echo "sv_bcbio=${bcbio_string}" >> $logfile