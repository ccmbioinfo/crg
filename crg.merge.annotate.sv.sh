#!/bin/bash
# usage: crg.merge.annotate.sv.sh <family>
# run from within family/bcbio-sv/ 

FAMILY=$1

#run from bcbio-sv folder
dir=`pwd`;
logfile="${dir}/${FAMILY}_sv_jobids.log";

if [ ! -f $logfile ]; then 
	touch $logfile;
fi;

#run metasv on each sample
for f in $FAMILY_*/$FAMILY/final/$FAMILY* 
do 
	cd $f
	SAMPLE="$(echo $f | cut -d'/' -f1)"
	metasv_jobs+=($(qsub ~/crg/metasv.pbs -v SAMPLE=$SAMPLE,FAMILY=$FAMILY))
	cd ../../../..
done
metasv_string=$( IFS=$':'; echo "${metasv_jobs[*]}" )

#run snpeff on each sample
for f in $FAMILY_*/$FAMILY/final/$FAMILY*
do
	SAMPLE="$(echo $f | cut -d'/' -f1)"
	snpeff_jobs+=($(qsub ~/crg/crg.snpeff.sh -F $f/$SAMPLE/variants.vcf.gz -W depend=afterany:"${metasv_string}"))
done
snpeff_string=$( IFS=$':'; echo "${snpeff_jobs[*]}" )


#run svscore on each sample 
for f in $FAMILY_*/$FAMILY/final/$FAMILY*
do
	SAMPLE="$(echo $f | cut -d'/' -f1)"
	svscore_jobs+=($(qsub ~/crg/crg.svscore.sh -F $f/$SAMPLE/variants.vcf.gz.snpeff.vcf -W depend=afterany:"${snpeff_string}"))
done
svscore_string=$( IFS=$':'; echo "${svscore_jobs[*]}" )

echo "metasv=${metasv_string}" >> ${logfile}
echo "snpeff=${snpeff_string}" >> ${logfile}
echo "svscore=${svscore_string}" >> ${logfile}
echo "Merging SVs with metaSV, annotating with snpeff, scoring with svscore, and creating report..." 

jobid=$(qsub ~/crg/crg.intersect_sv_vcfs.sh -F $FAMILY  -W depend=afterany:"${svscore_string}")

echo "intersect=$jobid">> ${logfile}



