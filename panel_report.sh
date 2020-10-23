#!/bin/bash

#usage: sh panel_report.sh <family> <path/to/ensemble.vcf.gz>
#run this from base <family> folder

function panel_report
{

if [ ! -d "${dir}/${family}" ]; then
	mkdir ${dir}/${family}
fi;

panel_vcf="${dir}/${family}/${family}-ensemble-annotated-decomposed.vcf.gz";
echo "$panel_vcf, $ensemble, $bed";

if [ ! -f "${panel_vcf}" ]; then
	echo "generating ${panel_vcf}";
	bedtools intersect -header -a $ensemble -b $bed > ${panel_vcf}
fi;

cd ${dir}
echo "report generation for $dir"
sh ~/cre/generate_reports.sh ${family} "wgs"
cd ..

}

family=$1;
ensemble=$2; #ensemble annotated vcf.gz from small-variant

#panel
dir="panel";
bed="genes/${family}.bed"
if [ ! -d genes ] || [ ! -f $bed ]; then
	echo "folder: genes or file: genes/${family}.bed not found. Exiting!";
	exit
fi;
panel_report

#panel flank
if [ ! -f "genes/${family}.flank.100k.bed" ]; then
	~/crg/crg.flank.sh $bed > genes/${family}.flank.100k.bed
fi;

dir="panel-flank100k";
bed="genes/${family}.flank.100k.bed";
panel_report
