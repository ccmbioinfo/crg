#PBS -N EHDN
#PBS -l vmem=80g,mem=80g,walltime=24:00:00
#PBS -joe
#PBS -d .


EHDN=/hpf/largeprojects/ccmbio/arun/Tools/EHDN.TCAG/ExpansionHunterDenovo-v0.7.0;
scripts=/hpf/largeprojects/ccmbio/arun/C4Rare/scripts/scripts;

if [ -z $family ]; then family=$1; fi;
outdir="str/expansion_hunter_denovo";

manifest="${outdir}/${family}_manifest.txt";
for i in `ls bcbio-align/${family}/final/${family}_*/${family}_*-ready.bam`; do

    if [ ! -d "$outdir" ]; then
        mkdir -p $outdir;
    fi;
    prefix=`basename $i -ready.bam`;s
    reads=`readlink -f $i`;
    ehdn_prefix="${outdir}/${prefix}";
    echo -e "/hpf/largeprojects/ccmbio/arun/Tools/EHDN.TCAG/ExpansionHunterDenovo-v0.7.0 --reference $ref --reads $reads --output-prefix $ehdn_prefix --min-anchor-mapq 50 --max-irr-mapq 40 \n"
    $EHDN --reference $ref --reads $reads --output-prefix $ehdn_prefix --min-anchor-mapq 50 --max-irr-mapq 40 &
    echo -e "${prefix}\tcase\t${ehdn_prefix}.json" >> $manifest;

done;

wait

##EHDN report per family
echo "generating multi-sample EHDN STR profile and report named ${outdir}/${family}_EHDN_str.tsv"
module load python/3.7.7

countsfile=${outdir}"/"${family}"_counts.txt";
output=${outdir}"/"${family}"_EHDN_str.tsv";
python ${scripts}/combine_counts.py --manifest $manifest --combinedCounts $countsfile
python ${scripts}/compare_anchored_irrs.py --manifest $manifest --inputCounts $countsfile --outputRegions $output --minCount 2 --testMet$
module unload python/3.7.7

