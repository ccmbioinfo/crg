#!/bin/bash
#PBS -l walltime=5:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g,mem=21g

vcf=$1

echo "Generating SV scores: " `date`
SVSCORE_DATA=/hpf/largeprojects/ccmbio/arun/Tools/SVScore
SVSCORE_SCRIPT=/hpf/largeprojects/ccmbio/naumenko/tools/svscore
module load perl/5.20.1

perl -w $SVSCORE_SCRIPT/svscore.pl -o max,sum,top5,top10,mean \
		 -e $SVSCORE_DATA/tests/refGene.exons.bed \
		 -f $SVSCORE_DATA/tests/refGene.introns.bed \
		 -dvc $SVSCORE_DATA/tests/whole_genome_SNVs.tsv.gz  \
		 -i $vcf > ${vcf}.svscore.vcf
