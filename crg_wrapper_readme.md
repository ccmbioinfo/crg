# CRG wrapper

This wrapper automates the crg pipeline starting from input preparation to running three main steps namely alignment (align), small variant calling (smv), and structural variant calling (sv) serially or each step independently if there is a need to re-run a specific step. It can process more than one project (family).

# Usage

crg_wrapper.py [-h] -f FILE -t {fastq,bam} -s {all,align,sv,smv} -d path

## Arguments

-h, --help show this help message and exit
-f FILE, --file FILE Three column TAB-seperated sample info file
-t {fastq,bam}, --type {fastq,bam}
File type to be used as input
-s {all,align,sv,smv}, --step {all,align,sv,smv}
Run a single step or all steps of crg pipeline
-d path, --dir path Absolute path to base directory where family directory
and bcbio sub-directories will be created

-f: Sample information for the project(s) is read from a three column TAB-seperated file. See ".tsv" files in `crg/crg_wrapper/test` folder.

- 1st column: This should be `<project>_<sampleid>`
- 2nd column: input path to CRAM file or read1 FASTQ(must have a "\_R1" pattern in filename). If there are more than one FASTQ file for read1, fastq column must be comma-seperated filenames for read 1 (currently, this expects that both the ends have same number of FASTQ files)
- 3rd column: input path to BAM file (cannot be more than one file)

When the sample info file has more than one project, the script assumes that the starting input file type ('-t bam' or '-t fastq') is same for all projects and expects the corresponding column to be non-empty. Following is an example sample info when fastq or bam column is non-empty.

| `sampleid` | `fastq`                                   | `bam` |
| ---------- | ----------------------------------------- | ----- |
| 111_AB     | 111_AB_R1.fastq.gz                        |       |
| 111_CD     | 111_CD_L001_R1.fq.gz,111_CD_L002_R2.fq.gz |       |
| 110_AF     | 110_AF_R1.fastq                           |       |
| 110_TF     | 110_TF_R1.fq                              |       |
| 200_AB     | 200_AB.cram                               |       |

| `sampleid` | `fastq` | `bam`      |
| ---------- | ------- | ---------- |
| 112_GF     |         | 112_GF.bam |

-t: The input type can be 'bam' or 'fastq' for align step and only 'bam' for sv and smv steps.

-s: All the options will use the input read from sample info file, and perform necessary directory creations/soft-linking using helper bash scripts. Reporting for smv and sv will be carried out only if a HPO file for that family is found in ~/gene_data/HPO folder.

- 'all': submits job for alignment first, followed by jobs for smv and sv in parallel
- 'smv': submits job for smv
- 'sv': submits job for sv

-d: Absolute path to directory where main `<family>` directory and other sub-directories will be created.

## Example

1. 'sv' and 'smv' can be run in parallel, since there are no conflict in directory structures

```bash
crg_wrapper.py -f 111.tsv -t bam -s smv -d /hpf/c4r/wgs
crg_wrapper.py -f 111.tsv -t bam -s sv -d /hpf/c4r/wgs
```

2. If running all three steps for a project, it is best to pass '-s all' instead of three seperate commands for each step as it takes care of job chaining (otherwise there will be error while align and sv is run in parallel, as they both depend on the `bcbio-align` directory)

```bash
crg_wrapper.py -f 111.tsv -t bam -s all -d /hpf/c4r/wgs
```

3. submit align/all jobs with 'fastq' input type

```bash
crg_wrapper.py -f 111.tsv -t fastq -s align -d /hpf/c4r/wgs
```

4. submit align/all jobs with 'bam' input type

```bash
crg_wrapper.py -f 111.tsv -t bam -s align -d /hpf/c4r/wgs
```

# To do

1. Use pandas to read sample info
2. If 'align' step for a project is already complete, provide option to run downstream steps from existing output rather than from sample info file
3. Add 'cram' to '-t' explicitly. Currently, it is read from fastq column and handled internally based on input extension to call `cram2fq.sh`
   - add calls to `cram2bam.sh`
   - reference is hard coded as /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa
4. Currently, not all steps of small-variant and structural variant reporting is fully automated (requires manual step). Integrate these steps when available.