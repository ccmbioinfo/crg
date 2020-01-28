import os
from os.path import abspath
import sys
import argparse
import uuid
import numpy
import pybedtools
import pysam
import subprocess

THREADS = 1
TMP_DIR=""
SAMPLE=""
METASV_OUTFILE="variants.vcf.gz"

def calc_paired_insert_stats(in_bam, nsample=1000000):
    """Taken straight from bcbio source code"""

    def insert_size_stats(dists):
        """Calcualtes mean/median and MAD from distances, avoiding outliers.
        MAD is the Median Absolute Deviation: http://en.wikipedia.org/wiki/Median_absolute_deviation
        """
        med = numpy.median(dists)
        filter_dists = list(filter(lambda x: x < med + 10 * med, dists))
        median = numpy.median(filter_dists)
        return {"mean": float(numpy.mean(filter_dists)), 
                "std": float(numpy.std(filter_dists)),
                "median": float(median),
                "mad": float(numpy.median([abs(x - median) for x in filter_dists]))}

    """Retrieve statistics for paired end read insert distances."""
    dists = []
    n = 0
    with pysam.Samfile(in_bam, "rb") as in_pysam:
        for read in in_pysam:
            if read.is_proper_pair and read.is_read1:
                n += 1
                dists.append(abs(read.isize))
                if n >= nsample:
                    break
    return insert_size_stats(dists)

def run_metasv(bam, manta, lumpy, wham):

    if not os.path.exists(SAMPLE):
        os.mkdir(SAMPLE)

    os.chdir(SAMPLE)

    if not os.path.exists(TMP_DIR):
        os.mkdir(TMP_DIR)

    insert_stats = calc_paired_insert_stats(bam)

    cmd = '''/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/run_metasv.py \
    --sample {} \
    --reference /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
    --bam {} \
    --outdir {} \
    --lumpy_vcf {} \
    --manta_vcf {} \
    --wham_vcf {} \
    --workdir {} \
    --num_threads {} \
    --assembly_max_tools=1 \
    --assembly_pad=500 \
    --isize_mean {} \
    --disable_assembly \
    --isize_sd {}'''.format(SAMPLE, bam, os.getcwd(), lumpy, manta, wham, TMP_DIR, THREADS, insert_stats['mean'], insert_stats['std'])

    subprocess.call(cmd, shell=True)

def filter_vcf():

    out_filtered = "{}.metasv.filtered.vcf.gz".format(SAMPLE)

    NEG_FILTER = ("((NUM_SVTOOLS = 1 && ABS(SVLEN)>50000) || "
            "(NUM_SVTOOLS = 1 && ABS(SVLEN)<4000 && BA_FLANK_PERCENT>80) || "
            "(NUM_SVTOOLS = 1 && ABS(SVLEN)<4000 && BA_NUM_GOOD_REC=0) || "
            "(ABS(SVLEN)<4000 && BA_NUM_GOOD_REC>2)) || FILTER='LowQual'")

    filter_cmd = '''bcftools filter \
    -e "{}" \
    -O z \
    -o {} \
    {}
    '''.format(NEG_FILTER, out_filtered, METASV_OUTFILE)

    subprocess.call(filter_cmd, shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs MetaSV to merge SV calls from various variant callers')
    parser.add_argument('-bam', type=str, help='Bam file of the sample', required=True)
    parser.add_argument('-manta', type=str, help='Manta vcf file')
    parser.add_argument('-lumpy', type=str, help='Lumpy vcf file')
    parser.add_argument('-wham', type=str, help='Wham vcf file')
    parser.add_argument('-sample', help='Sample file name', type=str, required=True)
    parser.add_argument('-threads', help='Number of threads, defaults to 4', type=int, default=4)
    args = parser.parse_args()

    if args.threads <= 0:
        ValueError("-threads needs to be at least 1")
    elif not args.sample:
        ValueError("-sample is empty")

    THREADS = args.threads

    TMP_DIR = "metasv-tmp-{}-{}".format(args.sample, uuid.uuid3(uuid.NAMESPACE_DNS, args.bam))

    SAMPLE = args.sample

    run_metasv(abspath(args.bam), abspath(args.manta), abspath(args.lumpy), abspath(args.wham))
    filter_vcf()