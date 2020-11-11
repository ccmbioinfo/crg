#!/bin/bash

#PBS -l walltime=4:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

# if [ -z $FAMILY ]; then
# 	echo "Specify family ID as first arguement to script"
# 	exit
# fi

# if [ -z $type ]; then
# 	echo "Specify type (filter or unfiltered) as second argument to script"
# 	exit
# fi

TODAY=`date +%Y-%m-%d`
FAMILY=$1

GENE_DATA=${HOME}/gene_data
HGMD=${GENE_DATA}/HGMD_2018/hgmd_pro.db
PROTEIN_CODING_GENES=${GENE_DATA}/grch37.p13.ensembl.sorted.protein.coding.genes.bed
EXON_BED=${GENE_DATA}/exons/hg19_UCSC_exons_canonical.bed
HPO=${GENE_DATA}/HPO/${FAMILY_ID}_HPO.txt
EXAC=${GENE_DATA}/ExAC/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt
OMIM=${GENE_DATA}/OMIM_2020-04-09/genemap2.txt
GNOMAD=${GENE_DATA}/gnomad_v2_sv.sites.bed
BIOMART=${GENE_DATA}/BioMaRt.GrCh37.75.ensembl.mim.hgnc.entrez.txt
MSSNG_MANTA_COUNTS=${GENE_DATA}/mssng_counts/Canadian_MSSNG_parent_SVs.Manta.counts.txt
MSSNG_LUMPY_COUNTS=${GENE_DATA}/mssng_counts/Canadian_MSSNG_parent_SVs.LUMPY.counts.txt

if [ ! $2 ]; then
	FAMILY=$1
	IN_FILES=`ls ${FAMILY}_*/${FAMILY}/final/${FAMILY}*/${FAMILY}*/*svscore.vcf | tr '\n' ' '`
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

for FILE in $IN_FILES
do
	bcftools filter -e "((NUM_SVTOOLS = 1 && ABS(SVLEN)>50000) || (NUM_SVTOOLS = 1 && ABS(SVLEN)<4000 && BA_FLANK_PERCENT>80) || (NUM_SVTOOLS = 1 && ABS(SVLEN)<4000 && BA_NUM_GOOD_REC=0) || (ABS(SVLEN)<4000 && BA_NUM_GOOD_REC>2)) || FILTER='LowQual'" \
		-O z \
		-o ${FILE}.filtered.gz \
		$FILE 
done

FILTERED_FILES=`ls ${FAMILY}_*/${FAMILY}/final/${FAMILY}*/${FAMILY}*/*filtered.gz | tr '\n' ' '`

for FILE in $FILTERED_FILES
do
	tabix $FILE
done

#make filtered report
echo "${PY} ${HOME}/crg/crg.intersect_sv_vcfs.py -protein_coding_genes=${PROTEIN_CODING_GENES} -exon_bed=${EXON_BED} -hgmd=${HGMD} -hpo=${HPO} -exac=${EXAC} -omim=${OMIM} -biomart=${BIOMART} -gnomad=${GNOMAD} -sv_counts ${MSSNG_MANTA_COUNTS} ${MSSNG_LUMPY_COUNTS} -o=${FAMILY}.wgs.sv.${TODAY}.tsv -i ${FILTERED_FILES}"
${PY} ${HOME}/crg/crg.intersect_sv_vcfs.py -protein_coding_genes=${PROTEIN_CODING_GENES} -exon_bed=${EXON_BED} -hgmd=${HGMD} -hpo=${HPO} -exac=${EXAC} -omim=${OMIM} -biomart=${BIOMART} -gnomad=${GNOMAD} -sv_counts ${MSSNG_MANTA_COUNTS} ${MSSNG_LUMPY_COUNTS} -o=${FAMILY}.wgs.sv.${TODAY}.tsv -i ${FILTERED_FILES}

${HOME}/crg/cap_report.py ${FAMILY}.wgs.sv.${TODAY}.tsv

if [[ "$OSTYPE" == *"darwin"* ]]; then
	open -a 'Microsoft Excel' ${FAMILY}.wgs.sv.${TODAY}.tsv
fi


#make unfiltered report
echo "${PY} ${HOME}/crg/crg.intersect_sv_vcfs.py -protein_coding_genes=${PROTEIN_CODING_GENES} -exon_bed=${EXON_BED} -hgmd=${HGMD} -hpo=${HPO} -exac=${EXAC} -omim=${OMIM} -biomart=${BIOMART} -gnomad=${GNOMAD} -sv_counts ${MSSNG_MANTA_COUNTS} ${MSSNG_LUMPY_COUNTS} -o=${FAMILY}.unfiltered.wgs.sv.${TODAY}.tsv -i ${IN_FILES}"
${PY} ${HOME}/crg/crg.intersect_sv_vcfs.py -protein_coding_genes=${PROTEIN_CODING_GENES} -exon_bed=${EXON_BED} -hgmd=${HGMD} -hpo=${HPO} -exac=${EXAC} -omim=${OMIM} -biomart=${BIOMART} -gnomad=${GNOMAD} -sv_counts ${MSSNG_MANTA_COUNTS} ${MSSNG_LUMPY_COUNTS} -o=${FAMILY}.unfiltered.wgs.sv.${TODAY}.tsv -i ${IN_FILES}

${HOME}/crg/cap_report.py ${FAMILY}.unfiltered.wgs.sv.${TODAY}.tsv

if [[ "$OSTYPE" == *"darwin"* ]]; then
	open -a 'Microsoft Excel' ${FAMILY}.unfiltered.wgs.sv.${TODAY}.tsv
fi
