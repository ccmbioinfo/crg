#!/bin/bash
# usage: crg.merge.annotate.sv.sh <family>
# run from within family/bcbio-sv/ 

FAMILY=$1

#run metasv on each sample
for f in $FAMILY_*/$FAMILY/final/$FAMILY* 
do 
	cd $f
	SAMPLE="$(echo $f | cut -d'/' -f1)"
	metasv_jobs+=($(qsub ~/crg/metasv.pbs -v SAMPLE=$SAMPLE,FAMILY=$FAMILY))
	cd ../../../..
done
metasv_string=$( IFS=$':'; echo "${metasv_jobs[*]}" )

#run snpeff on each sample on both filtered and unfiltered metasv output
#filtered metasv
for f in $FAMILY_*/$FAMILY/final/$FAMILY*
do
	SAMPLE="$(echo $f | cut -d'/' -f1)"
	snpeff_jobs+=($(qsub ~/crg/crg.snpeff.sh -F $f/$SAMPLE/*metasv.filtered.vcf.gz -W depend=afterany:"${metasv_string}"))
done
snpeff_string=$( IFS=$':'; echo "${snpeff_jobs[*]}" )
#unfiltered metasv
for f in $FAMILY_*/$FAMILY/final/$FAMILY*
do
	SAMPLE="$(echo $f | cut -d'/' -f1)"
	snpeff_unfiltered_jobs+=($(qsub ~/crg/crg.snpeff.sh -F $f/$SAMPLE/variants.vcf.gz -W depend=afterany:"${metasv_string}"))
done
snpeff_unfiltered_string=$( IFS=$':'; echo "${snpeff_unfiltered_jobs[*]}" )

#run svscore on each sample on both filtered and unfiltered metasv output
for f in $FAMILY_*/$FAMILY/final/$FAMILY*
do
	SAMPLE="$(echo $f | cut -d'/' -f1)"
	svscore_jobs+=($(qsub ~/crg/crg.svscore.sh -F $f/$SAMPLE/*filtered.vcf.gz.snpeff.vcf -W depend=afterany:"${snpeff_string}"))
done
svscore_string=$( IFS=$':'; echo "${svscore_jobs[*]}" )
#unfiltered metasv
for f in $FAMILY_*/$FAMILY/final/$FAMILY*
do
	SAMPLE="$(echo $f | cut -d'/' -f1)"
	svscore_unfiltered_jobs+=($(qsub ~/crg/crg.svscore.sh -F $f/$SAMPLE/variants.vcf.gz.snpeff.vcf -W depend=afterany:"${snpeff_unfiltered_string}"))
done
svscore_unfiltered_string=$( IFS=$':'; echo "${svscore_unfiltered_jobs[*]}" )
echo "Merging SVs with metaSV, annotating with snpeff, scoring with svscore, and creating report..."

qsub ~/crg/crg.intersect_sv_vcfs.sh -v family=$FAMILY,type=filtered -W depend=afterany:"${svscore_string}"
qsub ~/crg/crg.intersect_sv_vcfs.sh -v family=$FAMILY,type=unfiltered -W depend=afterany:"${svscore_unfiltered_string}"

