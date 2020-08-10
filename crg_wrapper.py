#!/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/python
# Python version >= 2.7
from __future__ import print_function
from collections import namedtuple
import sys, os, subprocess, glob, argparse, logging

# import pandas as pd

"""
ref links
---------
subprocess.check_output:
parse return value https://stackoverflow.com/questions/6657690/python-getoutput-equivalent-in-subprocess
buffer overflow https://stackoverflow.com/a/8700414
                https://docs.python.org/2.7/library/subprocess.html#subprocess.check_output

argparse:
custom argument type check https://stackoverflow.com/a/11541495
"""

# currently assumes location relative to user home
# move the bash scripts definition to a JSON file
user_home = os.path.expanduser("~")
ext_dict = {"fq": ".fq", "fastq": ".fq", "fastq.gz": ".fq.gz", "fq.gz": ".fq.gz"}
sampleinfo = namedtuple("sampleinfo", "familyid sampleid fastq bam")
crg_setup_script = os.path.expanduser("~/crg/crg.setup.sh")
bcbio_setup_script = os.path.expanduser("~/crg/crg.prepare_bcbio_run.sh")
bcbio_script = os.path.expanduser("~/cre/bcbio.pbs")
align_script = os.path.expanduser("~/crg/crg_wrapper/submit_align.sh")
smv_script = os.path.expanduser("~/crg/crg_wrapper/submit_smv.sh")
sv_script = os.path.expanduser("~/crg/crg_wrapper/submit_sv.sh")
cat_fastq = os.path.expanduser("~/crg/crg_wrapper/cat_fastq.sh")
cram2fq_script = os.path.expanduser("~/cre/cram2fq.sh")


def log_message(*message):

    """
    write message to logfile and stdout
    """

    if message:
        for i in message:
            logging.info(i)
            print(i)


def create_symlink(src, dest):

    if not os.path.isfile(dest):
        log_message("creating symlinks {} to {}".format(src, dest))
        subprocess.check_call(["ln", "-s", src, dest])


def concatenate_fastq(fastq):  # multiple fastq files per end

    """
    returns two list of fastq for concatenation
    """
    r1, r2 = [], []
    for read1 in fastq:
        prefix, suffix = read1.split("_R1")
        read2 = prefix + "_R2" + suffix
        if os.path.isfile(read1) and os.path.isfile(read2):
            r1.append(read1)
            r2.append(read2)
    return r1, r2


def check_fastq(fastq, sampleid, input_path):

    """
    based on input file suffix, performs
    1. conversion (cram2fq)
    2. concatenation (more than 1 fastq per end)
    3. create symlink if only one fastq per end
    """
    jobid = []
    if isinstance(fastq, list):  # multiple fastq files per end
        log_message("Multiple fastq files per end {}".format(sampleid))
        r1, r2 = concatenate_fastq(fastq)
        extension = r1[0].split("_R1")[1].split(".", 1)[1]
        qsub = ["qsub", cat_fastq, "-m", "ae", "-F"]
        x1 = (
            " ".join(r1)
            + " "
            + os.path.join(input_path, sampleid + "_1" + ext_dict[extension])
        )
        arg1 = ['"{}"'.format(x1)]
        x2 = (
            " ".join(r2)
            + " "
            + os.path.join(input_path, sampleid + "_2" + ext_dict[extension])
        )
        arg2 = ['"{}"'.format(x2)]
        qsub_cmd1 = qsub + arg1
        qsub_cmd2 = qsub + arg2
        log_message(qsub_cmd1, qsub_cmd2)
        cwd = os.getcwd()
        jobid1 = subprocess.check_output(qsub_cmd1).decode("UTF-8").rstrip()
        jobid2 = subprocess.check_output(qsub_cmd2).decode("UTF-8").rstrip()
        jobid = [jobid1, jobid2]

    elif ".cram" in fastq:
        log_message("CRAM input detected for {}".format(sampleid))
        cram = fastq
        if os.path.isfile(cram):
            cmd = "qsub {} -v cram={},sample={} -m ae".format(
                cram2fq_script, cram, os.path.join(input_path, sampleid)
            )
            log_message(cmd)
            jobid = [check_output(cmd.split(" ")).decode("UTF-8").rstrip()]
            log_message("cram2fq jobid for {}: {}".format(sampleid, jobid))

    elif ".fastq" in fastq or ".fq" in fastq:
        # fastq.split(".",1)[1] in ['fastq','fq','fastq.gz','fq.gz']: fails when path itself has "."
        log_message("Single FASTQ file per end for {} ".format(sampleid))
        prefix, suffix = fastq.split("_R1")
        extension = suffix.split(".", 1)[1]  # fq, fastq, fq.gz, fastq.gz
        r2 = prefix + "_R2" + suffix
        if os.path.isfile(fastq) and os.path.isfile(r2):
            dest = os.path.join(input_path, sampleid + "_1" + ext_dict[extension])
            create_symlink(fastq, dest)
            dest = os.path.join(input_path, sampleid + "_2" + ext_dict[extension])
            create_symlink(r2, dest)
    else:
        print(fastq)
        log_message(
            "Input {} given for {} is not handled. Exiting!".format(
                fastq.split(".", 1)[1], sampleid
            )
        )
        exit()
    return jobid


def prepare_input(bcbio_align, familyid, project, fq=1):

    """
    1. calls check_fastq to process fastq inputs
    2. creates symlink for bam in respective 'input' directories
    """
    input_path = os.path.join(bcbio_align, familyid, "input")
    jobid = []
    for samples in project:
        sampleid = samples.sampleid
        if fq == 1:
            if len(samples.fastq) > 0:
                fastq = samples.fastq
                log_message("Preparing input FASTQ files for {}".format(sampleid))
                j = check_fastq(fastq, sampleid, input_path)
                jobid += j
            else:
                log_message(
                    '"fastq" type given for {}, but fastq column is empty in sample info file. Exiting!'.format(
                        sampleid
                    )
                )
                exit()

        else:
            bam = samples.bam
            if not bam:
                log_message(
                    '"bam" type given for {}, but bam column is empty in sample info file. Exiting!'.format(
                        sampleid
                    )
                )
                exit()
            log_message(
                "Preparing/symlink input BAM files for {}".format(samples.familyid)
            )
            if os.path.isfile(bam):
                dest_link = os.path.join(input_path, sampleid + ".bam")
                create_symlink(bam, dest_link)
    return ":".join(jobid)


def setup_proj_directories(base_path, familyid):

    """
    creates directories required for crg pipeline using crg.setup.sh script
    """
    if not os.path.exists(os.path.join(base_path, familyid)):
        log_message("\n")
        log_message("Setting up crg directories for {}".format(familyid))
        cwd = os.getcwd()
        os.chdir(base_path)
        cmd = [crg_setup_script, familyid]
        subprocess.check_call(cmd)
        first = ["bcbio-align", "bcbio-sv", "bcbio-small_variants"]
        extra = ["genes", "panel", "panel-flank100k", "reports", "tcag"]
        os.chdir(cwd)
    else:
        log_message("\n")
        log_message(
            "crg directories for {}:{} already exists".format(
                familyid, os.path.join(base_path, familyid)
            )
        )


def submit_sv(
    base_path, align_path, sv_path, familyid, project, input_suffix, *prev_jobid
):

    """
    create neccessary directories for sv calling and submit jobs using submit_sv.sh script
    """
    if input_suffix == "bam":
        # align step was not done/not needed, start with bam file from sample_info
        if not all([os.path.exists(i) for i in [align_path, sv_path]]):
            setup_proj_directories(base_path, familyid)

        # adding this to ignore the 'bam' from sample_info file and start with bams found
        # in bcbio-align direcotry and create symlink only if bcbio-align final directory is empty
        # to do: must add --force argument to main script to make this check optional
        # glob won't raise error even if directory is missing
        bam_path = "{}_*/*.bam".format(
            os.path.join(align_path, familyid, "final", familyid)
        )
        if not glob.glob(bam_path):
            for i in project:
                bam_path = os.path.join(align_path, familyid, "final", i.sampleid)
                cmd = "mkdir -p {}".format(bam_path)
                # log_message(cmd)
                subprocess.check_call(cmd.split(" "))
                # os.makedirs(bam_path) #fails if some parents before leaf dir exists
                dest = os.path.join(bam_path, i.sampleid + "-ready.bam")
                create_symlink(i.bam, dest)
        else:
            log_message(
                "{} has final BAM files, skipping symlink step".format(align_path)
            )

    elif input_suffix == "":
        log_message(
            "structural variants calling will start after align job with jobid {} finishes".format(
                prev_jobid
            )
        )

    else:
        log_message(
            "input_suffix: {} not recognized, this must be bam if running '-s sv' from command-line or empty if part of '-s all'. Exiting!".format(
                input_suffix
            )
        )
        exit()

    if prev_jobid:  # first item in tuple is the jobid
        cmd = 'qsub {} -F "{}" -m ae -W depend=afterok:{}'.format(
            sv_script, familyid, prev_jobid[0]
        )
    else:
        cmd = 'qsub {} -F "{}" -m ae'.format(sv_script, familyid)

    cwd = os.getcwd()
    proj_dir = os.path.join(base_path, familyid)
    os.chdir(proj_dir)
    log_message("sv calling cmd:{}".format(cmd.split(" ")))
    jobid = subprocess.check_output(cmd.split(" ")).decode("UTF-8").rstrip()
    os.chdir(cwd)
    log_message("Job id for sv calling and report generation:{}\n".format(jobid))
    return jobid


def submit_smv(
    base_path, align_path, smv_path, familyid, project, input_suffix, *prev_jobid
):

    """
    create neccessary directories for smv calling and submit jobs using submit_smv.sh script
    """
    if input_suffix == "":
        # align step was already run and final bams are named with '-ready.bam'
        # must be renamed without '-ready' and symlinked to bcbio-small-variants
        log_message(
            "small-variant calling will start after align job with jobid {} finishes".format(
                prev_jobid[0]
            )
        )

    elif input_suffix == "bam":  # use the bam file from sample_info
        if not all(
            [os.path.exists(i) for i in [align_path, smv_path]]
        ):  # check not to overwrite here if align was run already
            log_message("project directories for crg not found. Creating..")
            setup_proj_directories(base_path, familyid)
        for i in project:
            dest = os.path.join(smv_path, familyid, "input", i.sampleid + ".bam")
            create_symlink(i.bam, dest)

    else:
        log_message(
            "input_suffix: {} not recognized, this must be bam if running '-s smv' from command-line or empty if part of '-s all'. Exiting!".format(
                input_suffix
            )
        )
        exit()

    if prev_jobid:
        cmd = 'qsub {} -F "{}" -W depend=afterok:{}'.format(
            smv_script, familyid, prev_jobid[0]
        )
    else:
        cmd = 'qsub {} -F "{}"'.format(smv_script, familyid)

    # submit jobs from base project/family directory
    cwd = os.getcwd()
    proj_dir = os.path.join(base_path, familyid)
    os.chdir(proj_dir)
    jobid = subprocess.check_output(cmd.split(" ")).decode("UTF-8").rstrip()
    os.chdir(cwd)
    log_message("smv calling cmd: {}".format(cmd.split(" ")))
    log_message("Job id for smv calling and report generation:{}\n".format(jobid))
    return jobid


def submit_align(base_path, align_path, familyid, project, input_suffix):

    """
    create neccessary directories for align and submit jobs using submit_sv.sh script
    """

    # setup script will create directories for crg pipeline and overwrite any
    # existing directories if found, so check for presence inside setup_proj_directories
    setup_proj_directories(base_path, familyid)
    cwd = os.getcwd()
    os.chdir(align_path)
    if input_suffix == "fastq":
        jobid = prepare_input(align_path, familyid, project, fq=1)
    elif input_suffix == "bam":
        jobid = prepare_input(align_path, familyid, project, fq=0)
    else:
        log_message("Unknown option for input suffix")
    os.chdir(cwd)

    if os.path.exists(align_path):
        cwd = os.getcwd()
        os.chdir(align_path)
        prep_qsub = "qsub {} -F".format(bcbio_setup_script)
        arg = ['"{} {}"'.format(familyid, "align_decoy")]
        # ~/crg/crg.prepare_bcbio_run.sh has no PBS directives and PBS_O_WORK = $HOME, so pass -d
        cmd = prep_qsub.split(" ") + arg + ["-d", align_path]
        if jobid:
            log_message("Job ids submitted for input prep: {}".format(jobid))
            cmd += ["-W depend=afterok:{}".format(jobid)]
        jobid = subprocess.check_output(cmd).decode("UTF-8").rstrip()
        log_message("bcbio setup cmd: {}".format(cmd.split(" ")))
        log_message("Job id for bcbio setup: {}".format(jobid))
        bcbio_qsub = [
            "qsub",
            bcbio_script,
            "-v",
            "project={}".format(familyid),
            "-m ae",
            "-W depend=afterok:{}".format(jobid),
        ]
        jobid = subprocess.check_output(bcbio_qsub).decode("UTF-8").rstrip()
        log_message("bcbio align cmd: {}".format(bcbio_qsub.split(" ")))
        log_message("Job id for bcbio align: {}\n".format(jobid))
        os.chdir(cwd)
    return jobid


def read_sample_info(filename):

    """
    reads three-cloumn tsv text file.
    todo: 
    - do validation check from args.file and 
    - read using pandas
    """

    projects = {}
    with open(filename) as f:
        # family_sampleid\tfastq\tbam
        column_names = f.readline().strip().split("\t")  # skip column names
        for i in f:
            sampleid, fastq, bam = i.strip("\n").split("\t")
            familyid = sampleid.split("_")[0]
            if "," in fastq:
                fastq = [items.strip() for items in fastq.split(",")]
            if not familyid in projects:
                projects[familyid] = []
            projects[familyid].append(sampleinfo(familyid, sampleid, fastq, bam))
    for i in projects:
        print(i, projects[i])
    return projects


def main(filename, input_suffix, step, base_path):

    projects = read_sample_info(filename)
    for familyid in projects:
        # setup project files crg.setup.sh familyid
        proj_dir = os.path.join(base_path, familyid)
        align_path = os.path.join(proj_dir, "bcbio-align")  # align with decoy
        sv_path = os.path.join(proj_dir, "bcbio-sv")  # structural variants
        smv_path = os.path.join(proj_dir, "bcbio-small-variants")  # small-variants
        project = projects[familyid]
        log_message("\nbase project dir for family {} : {}".format(familyid, proj_dir))
        if step == "all":
            # align
            align_jobid = submit_align(
                base_path, align_path, familyid, project, input_suffix
            )
            # input_suffix = ""
            # sv: structural variants
            sv_jobid = submit_sv(
                base_path, align_path, sv_path, familyid, project, "", align_jobid,
            )
            # smv: small variants
            smv_jobid = submit_smv(
                base_path, align_path, smv_path, familyid, project, "", align_jobid,
            )
        elif step == "align":
            # align
            jobid = submit_align(base_path, align_path, familyid, project, input_suffix)
        elif step == "sv":
            # sv: structural variants
            sv_jobid = submit_sv(
                base_path, align_path, sv_path, familyid, project, input_suffix
            )
        elif step == "smv":
            # smv: small variants
            smv_jobid = submit_smv(
                base_path, align_path, smv_path, familyid, project, input_suffix
            )
        else:
            log_message("Unknown option for step. Exiting!")
            exit()
    log_message("\nDONE!")


def valid_dir(dir):

    if os.path.isdir(dir):
        return dir
    else:
        message = "{} path does not exist. Please provide absolute path".format(dir)
        log_message(message)
        raise argparse.ArgumentTypeError(message)


def valid_file(filename):

    if not os.path.isfile(filename):
        message = "{} file does not exist".format(filename)
        log_message(message)
        raise argparse.ArgumentTypeError(message)
    else:
        if not os.path.getsize(filename) > 0:
            message = "{} file is empty".format(filename)
            log_message(message)
            raise argparse.ArgumentTypeError(message)
    return filename


if __name__ == "__main__":

    description = """This is a wrapper to prepare inputs and run crg steps: align, sv, smv
                    individually or all three based on a three column TAB-seperated sample_info file.
                    Example sample info files: "sample_info_wfq.tsv" and "sample_info_wbam.tsv" crg_wrapper/test directory.
                    Pipeline input is taken from 2nd (fastq) or 3rd (bam) column based on the option passed to -t or --type. 
                    Please make sure the corresponding column in the sample_info file is non-empty.
                    Detailed readme can be found here: crg/crg_wrapper_readme.md
                    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-f",
        "--file",
        type=valid_file,
        required=True,
        help="Three column TAB-seperated sample info file",
    )
    parser.add_argument(
        "-t",
        "--type",
        type=str,
        required=True,
        choices=["fastq", "bam"],
        help="File type to be used as input",
    )
    parser.add_argument(
        "-s",
        "--step",
        type=str,
        required=True,
        choices=["all", "align", "sv", "smv"],
        help="Run a single step or all steps of crg pipeline",
    )
    parser.add_argument(
        "-d",
        "--dir",
        type=valid_dir,
        required=True,
        metavar="path",
        help="Absolute path to base directory where family directory and bcbio sub-directories will be created",
    )
    args = parser.parse_args()
    logfile = args.file.split(".")[0] + "_" + args.step + ".LOG"
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="%(asctime)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M",
    )
    log_message("Sample info file: {}".format(args.file))
    log_message("File type to use as input: {}".format(args.type))
    log_message("Start crg pipeline from: {}".format(args.step))
    log_message("Create project directories under: {}".format(args.dir))
    log_message("Commands used are stored in log file: {}".format(logfile))
    main(args.file, args.type, args.step, args.dir)
