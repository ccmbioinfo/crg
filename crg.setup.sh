# script to set up a bcbio run

family=$1

if [ ! -d $family ]
then
	mkdir $family
	cd $family
	# create bcbio run directories
	bcbio_dirs=( "bcbio-align" "bcbio-small-variants" "bcbio-sv" )
	for dir in "${bcbio_dirs[@]}"
	do
		mkdir -p $dir/$family/input
	done
	# creat misc. dirs
	misc_dirs=( "genes" "panel" "panel-flank100k" "reports" "tcag" )
	for dir in "${misc_dirs[@]}"
	do
		mkdir $dir
	done
fi

echo "Report Directory Structure Set Up. Link the input files, named by familyID_sampleID.fq.gz and run next script"

