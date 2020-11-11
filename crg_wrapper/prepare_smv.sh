#!/bin/bash

##PBS -N prepare_snv
##PBS -joe
##PBS -l nodes=1:ppn=1
##PBS -l vmem=20g,mem=20g

#run this from base project/family directory
#usage: sh prepare_snv.sh 482

family_id=$1;
wd=`pwd`;

align_path="bcbio-align/${family_id}/final";
smv_path="bcbio-small-variants/${family_id}"
# echo "${align_path}";


smv_input=`ls ${smv_path}/input/*.bam | wc -l `;
echo "# smv input bams = ${smv_input}" #comment out
align_final=`ls ${align_path}/${family_id}_*/${family_id}*-ready.bam | wc -l`;
echo "# align final bams = ${align_final}"; #comment out


if [ "$smv_input" -eq "0" ]; then

	if [ "$align_final" -gt "0" ]; then
	
		for i in ${align_path}/${family_id}_*/${family_id}_*"-ready.bam"; do
			prefix=`basename $i -ready.bam`;
			src=`readlink -f $i`;
			dest="${smv_path}/input/${prefix}.bam";
			echo "creating symlink for $src -> $dest";
			ln -s $src $dest
		done;
	
	else
		echo "No BAM files found in ${align_path}. Exiting"
		exit
	fi;
#else part is symlink of BAM files are already present; proceed to dir setup
fi;


cd "$wd/bcbio-small-variants"
echo "Setting directories and config required for small vairants calling: ${family_id}";
~/crg/crg.prepare_bcbio_run.sh ${family_id} small_variants
cd "$wd"
	

