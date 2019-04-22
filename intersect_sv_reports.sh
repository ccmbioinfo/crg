#!/bin/bash

#PBS -l walltime=00:30:00,nodes=1:ppn=2
#PBS -joe .
#PBS -d .
#PBS -l vmem=5g,mem=5g

TODAY=`date +%Y-%m-%d`
if [ -z $1 ]; then
	echo "Specify family ID as first arguement to script"
	exit
fi
FAMILY_ID=$1
OUT=${FAMILY_ID}.wgs.sv.${TODAY}.tsv

GENE_DATA=${HOME}/gene_data
HGMD=${GENE_DATA}/HGMD_2018/hgmd_pro.db
EXON_BED=${GENE_DATA}/protein_coding_genes.exons.fixed.sorted.bed
HPO=${GENE_DATA}/HPO/${FAMILY_ID}_HPO.txt
EXAC=${GENE_DATA}/ExAC/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt
OMIM=${GENE_DATA}/OMIM_2018-11-01/genemap2.txt
GNOMAD=${GENE_DATA}/gnomad_v2_sv.sites.bed
BIOMART=${GENE_DATA}/BioMaRt.GrCh37.75.ensembl.mim.hgnc.entrez.txt

if [ ! $2 ]; then
	IN_FILES=`ls *.vcf | tr '\n' ' '`
else
	IN_FILES=${*:2}
fi

if [[ "$OSTYPE" == *"darwin"* ]]; then
	PY=python3
elif [[ "$OSTYPE" == "linux"* ]]; then
	module load python/3.7.1
	module load bedtools
	module load sqlite/3.20.0

	PY=python
fi

echo "${PY} ${HOME}/crg/intersect_sv_vcfs.py -exon_bed=${EXON_BED} -hgmd=${HGMD} -hpo=${HPO} -exac=${EXAC} -omim=${OMIM} -biomart="${BIOMART} -gnomad=${GNOMAD}" -o=${OUT} -i ${IN_FILES}"
${PY} ${HOME}/crg/intersect_sv_vcfs.py -exon_bed=${EXON_BED} -hgmd=${HGMD} -hpo=${HPO} -exac=${EXAC} -omim=${OMIM} -biomart=${BIOMART} -gnomad=${GNOMAD} -o=${OUT} -i ${IN_FILES}

if [[ "$OSTYPE" == *"darwin"* ]]; then
	open -a 'Microsoft Excel' ${OUT}
fi
