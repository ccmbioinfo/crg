#!/bin/bash

#PBS -N submit_smv
#PBS -joe
#PBS -d .
#PBS -l nodes=1:ppn=1
#PBS -l vmem=20g,mem=20g
#PBS -m ae

#submit this from base project directory: /hpf/largeprojects/ccmbio/aarthi/proj_CHEO/CRG/496
#usage: qsub submit_smv.sh -F "496"


family_id=$1;
wd=`pwd`;
logfile="${wd}/bcbio-small-variants/${family_id}_smv_jobids.log"; #all job ids are recorded here
if [ -f $logfile ]; then	
	rm $logfile
fi;
touch ${logfile};

#step1: symlink input files and setup bcbio directories
echo "curr dir from submit_smv: $wd";
if [ -d "bcbio-align" ] && [ -d "bcbio-small-variants" ]; then 
    sh ~/crg/crg_wrapper/prepare_smv.sh ${family_id}
else
    echo "bcbio-small-variants and bcbio-align directories not found."
    echo "Run ~/crg/crg.setup.sh to setup main directories before running this"
    exit
fi


#step2: submit bcbio jobs for small-variant calling
#do the following steps from within bcbio-small-variants directory
cd bcbio-small-variants
bcbio_job=$(qsub ~/cre/bcbio.pbs -v project=${family_id});
echo "submitted smv bcbio jobid for ${family_id}: ${bcbio_job}";

#step3: cleanup (got walltime(30H) exceeeded error for one sample, need to check this)
cleanup_job=$(qsub ~/cre/cre.sh -v family=${family_id},cleanup=1,make_report=0,type=wgs -W depend=afterok:${bcbio_job});
echo "submitted cleanup jobid for ${family_id}: ${cleanup_job}";

#step4: small variant reporting
smv_report=$(qsub ~/cre/cre.sh -v family=${family_id} -W depend=afterok:${cleanup_job});
echo "submitted smv reporting jobid for ${family_id}: ${smv_report}";

#record submitted jobids here
echo "smv_bcbio=${bcbio_job}" >> ${logfile};
echo "smv_cleanup=${cleanup_job}" >> ${logfile};
echo "smv_report=${smv_report}" >> ${logfile};

#panel steps: needs the panel and panel100k to be ready. 
#for now run it seperately 


