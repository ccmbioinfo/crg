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

#get manta BNDs supported by at least 5 split reads or read pairs
for f in $FAMILY_*/$FAMILY/final/$FAMILY*
do
	bcftools view -i 'INFO/SVTYPE="BND" && (FORMAT/PR[0:1] >= 5 || FORMAT/SR[0:1] >= 5)' -O z $f/${FAMILY}-manta.vcf.gz > $f/${FAMILY}-manta.BND.vcf.gz
done


#run snpeff on manta vcf for each sample
for f in $FAMILY_*/$FAMILY/final/$FAMILY*
do
	snpeff_manta_jobs+=($(qsub ~/crg/crg.snpeff.sh -F $f/${FAMILY}-manta.BND.vcf.gz))
done
snpeff_manta_string=$( IFS=$':'; echo "${snpeff_manta_jobs[*]}" )


#run svscore on manta vcf for each sample 
for f in $FAMILY_*/$FAMILY/final/$FAMILY*
do
	SAMPLE="$(echo $f | cut -d'/' -f1)"
	svscore_manta_jobs+=($(qsub ~/crg/crg.svscore.sh -F $f/${FAMILY}-manta.BND.vcf.gz.snpeff.vcf -W depend=afterany:"${snpeff_manta_string}"))
done
svscore_manta_string=$( IFS=$':'; echo "${svscore_manta_jobs[*]}" )


#run metasv on each sample
for f in $FAMILY_*/$FAMILY/final/$FAMILY* 
do 
	cd $f
	SAMPLE="$(echo $f | cut -d'/' -f1)"
	metasv_jobs+=($(qsub ~/crg/metasv.pbs -v SAMPLE=$SAMPLE,FAMILY=$FAMILY))
	cd ../../../..
done
metasv_string=$( IFS=$':'; echo "${metasv_jobs[*]}" )

#run snpeff on metasv vcf for each sample
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

depend_string=`echo ${svscore_string}:${svscore_manta_string}`
jobid=$(qsub ~/crg/crg.intersect_sv_vcfs.sh -F $FAMILY  -W depend=afterany:"${depend_string}")

echo "intersect=$jobid">> ${logfile}



