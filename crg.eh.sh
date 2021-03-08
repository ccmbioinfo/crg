#PBS -N EH
#PBS -l vmem=80g,mem=80g,walltime=24:00:00
#PBS -joe
#PBS -d .

ref=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37d5/seq/GRCh37d5.fa;
catalog=/hpf/largeprojects/ccmbio_ephemeral/ExpansionHunter/tandem_repeat_disease_loci_v1.1.hg19.masked.json;
scripts=/hpf/largeprojects/ccmbio/aarthi/proj_CHEO/CRG/str;

#original annotation file from Brett:/hpf/largeprojects/ccmbio_ephemeral/ExpansionHunter/tandem_repeat_disease_loci_v1.1.tsv
#em_dashes from above tsv were causing issues when splitting annotation in eh_sample_report.py; so manually replaced those with hyphen and stored in dir below
tr_annot=/hpf/largeprojects/ccmbio/aarthi/proj_CHEO/CRG/str/tandem_repeat_disease_loci_v1.1.tsv;


if [ -z $family ]; then family=$1; fi;

bcbio="${family}/bcbio-align/${family}/final/${family}_*";
if [ -z "`ls ${bcbio}/${family}_*-ready.bam`" ]; then 
    if [ -z "`ls ${family}/${family}*-ready.bam`" ]; then
        echo "BAM files for $family not found inside bcbio-align/ and $family/. Exiting!";
        exit
    else
        dir=$family;
    fi;
else
    dir=$bcbio;
fi;

outdir="${family}/str/expansion_hunter";
echo "${outdir}, ${dir}, $family";

module load ExpansionHunter/3.0.1
for i in `ls $dir/${family}*-ready.bam`; do
    echo "bam = $i";
    if [ ! -d $outdir ]; then
        mkdir -p $outdir;
    fi;
    prefix=`basename $i -ready.bam`;
    reads=`readlink -f $i`;
    eh_prefix="${outdir}/${prefix}";
    sex=`sh ~/crg/get_XY.sh $reads`;
    echo  -e "COMMAND: ExpansionHunter --reads $reads --reference $ref --variant-catalog $catalog --output-prefix $eh_prefix --sex $sex \n"
    ExpansionHunter --reads $reads --reference $ref --variant-catalog $catalog --output-prefix $eh_prefix --sex $sex 
done
module unload ExpansionHunter/3.0.1

##EH report per family
echo "generating multi-sample EH STR report named ${outdir}/${family}_EH_str.tsv"

#following Python script requires csv, python-docx
#modules I checked didn't have python-docx, so installed in my conda env "str"
source /hpf/largeprojects/ccmbio/aarthi/miniconda3/etc/profile.d/conda.sh
conda activate str
tsv="${outdir}/${family}_EH_str.tsv";
echo "COMMAND: python ${scripts}/generate_EH_genotype_table.generic.py ${outdir}"
python ${scripts}/generate_EH_genotype_table.generic.py ${outdir} > $tsv

#annotate with ref gene, location and repeat size
annot="${outdir}/${family}_EH_str.annot.tsv"
echo "annotating multi-sample EH STR report named $annot"
echo "COMMAND: python ${scripts}/add_gene+threshold_to_EH_column_headings2.py $tsv ${tr_annot} > $annot"
python ${scripts}/add_gene+threshold_to_EH_column_headings2.py $tsv ${tr_annot} > $annot

#transpose and seperate columns to final report
xlsx="${outdir}/${family}_EH_v1.1.xlsx";
echo "generating family-wise report: $xlsx"
echo "COMMAND: python ~/crg/eh_sample_report.py ${annot} ${xlsx}"
python ~/crg/eh_sample_report.py ${annot} ${xlsx}



