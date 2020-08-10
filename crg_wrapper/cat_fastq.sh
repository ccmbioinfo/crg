#!/bin/bash

#PBS -N concat_fastq
#PBS -d .
#PBS -joe
#PBS -l vmem=20g,mem=20g
#PBS -l walltime=10:00:00,nodes=1:ppn=1

#usage: sh cat_fastq.sh file1 file2 file3 outfile
#takes any number of input files for concatenation 
#last argument is treated as outputfile
#program quits if either any of the input file is not found or output file already exists

#shift: https://stackoverflow.com/questions/32681930/bash-shell-passing-variable-number-of-arguments
#$@ -> arguments passed to script (bash fails to split arguments by space when called from Python, so shift command won't work; use read command ) 


inputfiles=();
inputfiles=(`echo $@ | awk '{ for(i=1;i<NF;i++) print $i; }'`);
outfile=`echo $@ | awk '{ print $NF; }'`; #last arg is output file

for i in ${inputfiles[@]}; do
	if [ ! -f $i ]; then 
		echo "input file: $i not found. Exiting!";
		exit
	fi;
done;

if [ -f "$outfile" ]; then
	echo "output file: $outfile already exists, not running concatenation. Exiting!";
	exit
fi;

echo "inputfiles = ${inputfiles[*]}"
echo "outfile = $outfile"

cat_string=$(IFS=" " ; echo "${inputfiles[*]}");
echo "Running: cat ${cat_string} > ${outfile}";
cat ${cat_string} > ${outfile};
echo "Done";