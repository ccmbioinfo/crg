#!/bin/bash


#readme: use this in the main submission scripts after bcbio-align step has run
#usage: run this from within <family> directory
#sh crg_str.sh <family> 

family=$1;

#ppn=`ls bcbio-align/${family}/final/${family}_*/${family}_*-ready.bam | wc -l`;

#detect repeats in known sites in each sample, combine and create a annotated xlsx report per family
eh_job=$(qsub ~/crg/crg.eh.sh -v family=$family);
echo "Submitted EH for $family: ${eh_job}";

#detect denovo repeats in each sample, combine and create report per family
ehdn_job=$(qsub ~/crg/crg.ehdn.sh -v family=$family,pipeline=crg);
echo "Submitted job EHDN: ${ehdn_job}"

outdir="${family}/str/expansion_hunter_denovo"; #for bcbio-runs, change this if using this script elsewhere
qsub ~/crg/ehdn_report.sh -F "${family} ${outdir}" -W depend=afterok:${ehdn_job}
echo "Submitted job for DBSCAN clustering and report"


