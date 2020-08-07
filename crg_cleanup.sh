#!/bin/bash
#this script is for cleaning up a completed crg run prior to moving
#results to /hpf/largeprojects/ccm_dccforge/dccdipg/c4r_wgs/results/
#usage: crg_cleanup.sh <family>
#run from within top level of crg directory

FAMILY=$1

#move reports into reports directory
mv bcbio-small-variants/${FAMILY}/${FAMILY}.wes.regular*csv reports
mv bcbio-sv/${FAMILY}.wgs.sv*tsv reports
mv tcag/${FAMILY}*cnv.tsv reports

#rename panel report and move into reports directory
PANEL="$(basename panel/${FAMILY}/${FAMILY}.wgs*csv)"
PANEL_DATE="$(echo $PANEL | cut -d'.' -f3)"

mv panel/${FAMILY}/${FAMILY}.wgs*csv panel/${FAMILY}/${FAMILY}.wgs.panel.${PANEL_DATE}.csv
mv panel/${FAMILY}/${FAMILY}.wgs.panel.${PANEL_DATE}.csv reports

#rename panel flank report and move into reports directory
PANEL_FLANK="$(basename panel-flank100k/${FAMILY}/${FAMILY}.wgs*csv)"
PANEL_FLANK_DATE="$(echo $PANEL_FLANK | cut -d'.' -f3)"

mv panel-flank100k/${FAMILY}/${FAMILY}.wgs*csv  panel-flank100k/${FAMILY}/${FAMILY}.wgs.panel-flank100k.${PANEL_FLANK_DATE}.csv
mv panel-flank100k/${FAMILY}/${FAMILY}.wgs.panel-flank100k.${PANEL_FLANK_DATE}.csv reports

#move bams and multiqc to top level directory
mv bcbio-align/${FAMILY}/final/${FAMILY}_*/*bam* .
mv bcbio-small-variants/${FAMILY}/multiqc .

#remove unecessary folders
rm -r  bcbio-align/${FAMILY}/final/  bcbio-align/${FAMILY}/work/
rm -r  bcbio-small-variants/${FAMILY}/sv/
rm -r remove_decoys

#cleanup bcbio-sv directory
rm bcbio-sv/${FAMILY}*/${FAMILY}/final/${FAMILY}_*/*bam*
rm -r bcbio-sv/${FAMILY}*/${FAMILY}/work
rm -r bcbio-sv/${FAMILY}*/${FAMILY}/final/${FAMILY}_*/qc
rm -r bcbio-sv/svscoretmp bcbio-sv/*crg* bcbio-sv/snpEff_genes.txt bcbio-sv/snpEff_summary.html bcbio-sv/conf.toml

#cleanup genes directory
rm genes/*missing* genes/*unsorted*
