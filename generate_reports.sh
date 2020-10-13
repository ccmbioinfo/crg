#!/bin/bash
# generate both the regular and synonymous reports for a family
# pass in the family id as the first positional parameter, use -email to recieve an email when both jobs are done
# usage: generate_reports.sh <family id> <report_type> [optional -email] 

family=$1
report_type=$2


if [ "$3" == "-email" ]; then
  email_flag="-m e"
else
	email_flag=""
fi

cd $family
family_vcf="${family}-ensemble-annotated-decomposed.vcf.gz"

if [ -f $family_vcf ]; then
	if [ "$report_type" = wes ]; then
		vcf2cre_job="$(qsub ~/cre/cre.vcf2cre.sh -v original_vcf="${family}-gatk-haplotype-annotated-decomposed.vcf.gz",project=${family})"
	elif [ "$report_type" = wgs ]; then
		vcf2cre_job="$(qsub ~/crg/crg.vcf2cre.sh -v original_vcf="${family}-ensemble-annotated-decomposed.vcf.gz",project=${family})"
	fi
else
	echo "${family_vcf} Not Present, Aborting."
	cd ..
	exit
fi

if [ "$report_type" = wes ]; then
	standard_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${vcf2cre_job}" -v family=${family})"
	echo "Standard WES Report Job ID: ${standard_job}"
elif [ "$report_type" = wgs ]; then
	standard_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${vcf2cre_job}" -v family=${family},type=wgs)"
        echo "WGS Report Job ID: ${standard_job}"
fi


echo "The Rerun subfolder will be renamed by the current date after the reports are created"
cleanup_job="$(qsub ~/cre/rename_rerun.sh -W depend=afterok:"${standard_job}" -v family=${family})"

cd ..
