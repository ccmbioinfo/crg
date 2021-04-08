#PBS -N EHDN_REPORT
#PBS -l vmem=80g,mem=80g,walltime=48:00:00
#PBS -joe
#PBS -d .


scripts=~/crg/str;
EHDN_files=/hpf/largeprojects/ccm_dccforge/dccdipg/Common/annotation/ExpansionHunterDenovo;
g1k_outlier="${EHDN_files}/1000G_outlier";
omim=${EHDN_files}/OMIM_hgnc_join_omim_phenos_2020-08-05.tsv;
gnomad=${EHDN_files}/gnomad.v2.1.1.lof_metrics.by_gene.txt;
annovar=/hpf/largeprojects/ccm_dccforge/dccdipg/Common/crg2-non-conda-tools/annovar;
annovar_db=/hpf/largeprojects/ccm_dccforge/dccdipg/Common/annotation/annovar/humandb/;

family=$1;
outdir=$2;

##find outlier(noise) samples in current cohort to be removed later
echo "STEP1: write outlier samples to a file"
profile="${outdir}/${family}_EHDN_str.tsv";
outliers="${outdir}/${family}_outliers.txt";
python ${scripts}/find_outliers.py ${profile} ${outliers}
cat ${outliers} ${g1k_outlier} > temp && mv temp ${outliers}

##run outlier detection using DBSCAN clustering with 1000G as control
echo "STEP2: Outlier detection using DBSCAN clustering"
outpath="${outdir}/outliers"
Rscript ${scripts}/DBSCAN.EHdn.parallel.R --infile ${profile} --outpath ${outpath} --outlierlist ${outliers}
outlier_tsv=`ls ${outpath}/EHdn.expansions*tsv`;
Rscript ${scripts}/mergeExpansions.R --ehdn ${profile}  --outlier ${outlier_tsv} --outpath ${outpath}


#following Python script requires python3 and xlsxwriter
#installed in my conda env "str"
source /hpf/largeprojects/ccmbio/aarthi/miniconda3/etc/profile.d/conda.sh
conda activate str
##format above output for annovar annotation
echo "STEP3.1: format DBSCAN output for ANNOVAR annotation"
merge_tsv=`ls ${outpath}/merged.rare.expansions.[0-9]*.tsv`;
python ${scripts}/format_for_annovar.py ${outlier_tsv} ${merge_tsv}

echo "STEP3.2: annotate with OMIM"
annovar_input=`ls ${outpath}/merged.rare.EHdn.expansion.[0-9]*tsv`;
outfile=`echo ${merge_tsv} | awk '{ split($1,a,".tsv"); print a[1]; }'`;
${annovar}/table_annovar.pl ${annovar_input} ${annovar_db} -buildver hg19 -outfile "${outfile}-OMIM" -remove --onetranscript --otherinfo -protocol refGene -operation gx -nastring . -xreffile ${omim}

echo "STEP3.3: annotate with gnoMAD"
${annovar}/table_annovar.pl ${annovar_input} ${annovar_db} -buildver hg19 -outfile "${outfile}-gnoMAD" -remove --onetranscript --otherinfo -protocol refGene -operation gx -nastring . -xreffile ${gnomad}

echo "STEP3.4: format annotated reports to xlsx: EHDN.rare.merged.annotated.xlsx"
gnomad_out=`ls ${outfile}-gnoMAD.hg19_multianno.txt`
omim_out=`ls ${outfile}-OMIM.hg19_multianno.txt`
date=`date +%Y-%m-%d`
xlsx="${outpath}/${family}.EHDN.${date}.xlsx"
python ${scripts}/format_from_annovar.py ${gnomad_out} ${omim_out} ${xlsx}