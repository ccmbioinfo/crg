#!/bin/bash

# input is a family id and the location of the previous run
# run this script in a place with enough space and within an empty folder

family_id=$1
project_folder=$2

if [ ! -d ${project_folder} ]
then
	echo "Family folder: " $project_folder "Not Found"
	exit
fi

if [ -d "remove_decoys" ]
then
	echo "Decoys folder already exists."
	exit
fi

if [ -d "sv-calls" ]
then
	echo "sv-calls folder already exists"
	exit
fi

echo "Linking bam files from ${project_folder}"
mkdir remove_decoys
cd remove_decoys
ln -s ${project_folder}/*.bam .
cd ..

echo "Creating directory structure for individual SV calling"
mkdir sv-calls
cd sv-calls
for f in ../remove_decoys/*.bam; 
do
	mkdir -p $(echo $(basename $f)/"${family_id}"/input | sed 's/\.bam//g'); 
done

cd ../remove_decoys

echo "Submitting Jobs to Remove Decoy Reads"
decoy_jobs=()
for f in *.bam;
do
	decoy_jobs+=($(qsub ~/cre/cre.bam.remove_decoy_reads.sh -v bam=$f))
done
decoy_string=$( IFS=$':'; echo "${decoy_jobs[*]}" ) # store job ids

echo "Preparing bcbio SV calling (will run once decoys removed)"
cd ../sv-calls
setup_bcbio=$(qsub ~/crg/crg.setup.individual-sv.sh -v family_id="${family_id}" -W depend=afterany:"${decoy_string}")

echo "Submitting bcbio runs (will start once runs are prepared)"
bcbio_jobs=()
for f in *; 
do 
	cd $f;
	bcbio_jobs+=($(qsub ~/cre/bcbio.pbs -v project="${family_id}" -W depend=afterany:"${setup_bcbio}"))
	cd ..; 
done
bcbio_string=$( IFS=$':'; echo "${bcbio_jobs[*]}" )

echo "Setup Complete."
echo "Once bcbio jobs: " $( IFS=$', '; echo "${bcbio_jobs[*]}" ) " are done, run the following commands from within the sv-calls folder"
echo "MetaSV:"
echo '	for f in '${family_id}'_*/'${family_id}'/final/'${family_id}'*; do cd $f; pwd; ls; qsub ~/crg/metasv.pbs -v PROJECT=${family_id},SAMPLE="$(echo $f | sed -n -e 's/.*_//p')"; cd ../../../..; done")'
echo "SnpEff:"
echo '	for f in '${family_id}'_*/'${family_id}'/final/'${family_id}'_*/*/*metasv*.gz; do qsub ~/crg/crg.snpeff.sh -F $f; ls; done'
echo "SVScore:"
echo '	for f in '${family_id}'_*/'${family_id}'/final/'${family_id}'_*/*/*snpeff*.vcf; do qsub ~/crg/crg.svscore.sh -F $f; done'
echo "Then create a directory to combine SVScore VCFs and link them:"
echo '	mkdir combine_vcfs'
echo '	cd combine_vcfs'
echo '	for f in ../'${family_id}'_*/'${family_id}'/final/'${family_id}'_*/*/*svscore*.vcf;do ln -s $f . ; done'
echo "Finally, run the intersect sv script:"
echo '	~/crg/crg.intersect_sv_vcfs.sh -F '${family_id}