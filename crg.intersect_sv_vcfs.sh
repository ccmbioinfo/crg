#!/bin/bash

#PBS -l walltime=0:30:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $1 ]; then
	echo "Specify family ID as first arguement to script"
	exit
fi

TODAY=`date +%Y-%m-%d`
FAMILY_ID=$1
OUT=${FAMILY_ID}.wgs.sv.${TODAY}.tsv

GENE_DATA=${HOME}/gene_data
HGMD=${GENE_DATA}/HGMD_2018/hgmd_pro.db
PROTEIN_CODING_GENES=${GENE_DATA}/grch37.p13.ensembl.sorted.protein.coding.genes.bed
EXON_BED=${GENE_DATA}/protein_coding_genes.exons.fixed.sorted.bed
HPO=${GENE_DATA}/HPO/${FAMILY_ID}_HPO.txt
EXAC=${GENE_DATA}/ExAC/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt
OMIM=${GENE_DATA}/OMIM_2018-11-01/genemap2.txt
GNOMAD=${GENE_DATA}/gnomad_v2_sv.sites.bed
BIOMART=${GENE_DATA}/BioMaRt.GrCh37.75.ensembl.mim.hgnc.entrez.txt
MSSNG_MANTA_COUNTS=${GENE_DATA}/mssng_counts/Canadian_MSSNG_parent_SVs.Manta.counts.txt
MSSNG_LUMPY_COUNTS=${GENE_DATA}/mssng_counts/Canadian_MSSNG_parent_SVs.LUMPY.counts.txt

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

echo "${PY} ${HOME}/crg/crg.intersect_sv_vcfs.py -protein_coding_genes=${PROTEIN_CODING_GENES} -exon_bed=${EXON_BED} -hgmd=${HGMD} -hpo=${HPO} -exac=${EXAC} -omim=${OMIM} -biomart=${BIOMART} -gnomad=${GNOMAD} -sv_counts ${MSSNG_MANTA_COUNTS} ${MSSNG_LUMPY_COUNTS} -o=${OUT} -i ${IN_FILES}"
${PY} ${HOME}/crg/crg.intersect_sv_vcfs.py -protein_coding_genes=${PROTEIN_CODING_GENES} -exon_bed=${EXON_BED} -hgmd=${HGMD} -hpo=${HPO} -exac=${EXAC} -omim=${OMIM} -biomart=${BIOMART} -gnomad=${GNOMAD} -sv_counts ${MSSNG_MANTA_COUNTS} ${MSSNG_LUMPY_COUNTS} -o=${OUT} -i ${IN_FILES}

${HOME}/crg/cap_report.py ${OUT}

if [[ "$OSTYPE" == *"darwin"* ]]; then
	open -a 'Microsoft Excel' ${OUT}
fi
