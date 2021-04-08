#PBS -N EH
#PBS -l vmem=80g,mem=80g,walltime=24:00:00
#PBS -joe
#PBS -d .

ref=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37d5/seq/GRCh37d5.fa;
catalog=/hpf/largeprojects/ccmbio/arun/C4Rare/Tandem_repeat_disease_loci.hg19.json;
scripts=/hpf/largeprojects/ccmbio/arun/C4Rare/scripts/scripts;
tr_annot=/hpf/largeprojects/ccmbio/arun/C4Rare/scripts/scripts/TR_disease_loci_for_adding_to_col_headings.tsv;


if [ -z $family ]; then family=$1; fi;
outdir="str/expansion_hunter";

module load ExpansionHunter/3.0.1
for i in `ls bcbio-align/${family}/final/${family}_*/${family}_*-ready.bam`; do
    if [ ! -d $outdir ]; then
        mkdir -p $outdir;
    fi;
    prefix=`basename $i -ready.bam`;
    reads=`readlink -f $i`;
    eh_prefix="${outdir}/${prefix}";
    sex=`sh ~/crg/get_XY.sh $reads`;
    echo  -e "ExpansionHunter --reads $reads --reference $ref --variant-catalog $catalog --output-prefix $eh_prefix --sex $sex \n"
    ExpansionHunter --reads $reads --reference $ref --variant-catalog $catalog --output-prefix $eh_prefix --sex $sex &
done
module unload ExpansionHunter/3.0.1

##EH report per family
echo "generating multi-sample EH STR report named ${outdir}/${family}_EH_str.tsv"

#following Python scripts requires csv, python-docx
#modules I checked didn't have python-docx, so installed in my conda env "str"
source /hpf/largeprojects/ccmbio/aarthi/miniconda3/etc/profile.d/conda.sh
conda activate str
tsv="${outdir}/${family}_EH_str.tsv";
python ${scripts}/generate_EH_genotype_table.generic.py ${eh_prefix} > $tsv

#annotate with ref gene, location and repeat size
python ${scripts}/add_gene+threshold_to_EH_column_headings2.py $tsv ${tr_annot}
annot="${outdir}/${family}_EH_str.annot.tsv"
xlsx="${outdir}/${family}_EH_str.annot.xlsx"
#repeat GT in tsv gets auto-formatted to date in Excel, so convert tsv to xlsx
python ${scripts}/tsv2xlsx.py $annot $xlsx
