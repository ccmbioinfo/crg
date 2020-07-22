import argparse
import numpy as np
import pandas as pd
from pybedtools import BedTool


def intersect_variants(cnv, sv):
    # get coordinates of DEL and DUP from SV and CNV reports
    sv_coord = sv[['CHROM', 'POS', 'END', 'SVTYPE']]
    sv_coord = sv_coord[(sv_coord['SVTYPE'] != 'INV')
                        & (sv_coord['SVTYPE'] != 'INS')]
    cnv_coord = cnv[['CHROM', 'START', 'END', 'SVTYPE']]

    # convert report dataframes to bedtools
    cnv_bed = BedTool.from_dataframe(cnv_coord)
    sv_bed = BedTool.from_dataframe(sv_coord)

    # intersect CNV and SV bed files with reciprocal overlap of 50%
    intersection = cnv_bed.intersect(sv_bed, wa=True, wb=True, F=0.5,
                                     f=0.5).saveas('intersection.bed')
    intersection = pd.read_csv('intersection.bed', sep='\t', names=['CNV_CHROM', 'CNV_START', 'CNV_END', 'CNV_SVTYPE', 'SV_CHROM',
                                                                    'SV_START', 'SV_END', 'SV_SVTYPE'])
    # remove overlapping SVs that are different SVTYPES
    intersection = intersection[intersection['CNV_SVTYPE'] == intersection['SV_SVTYPE']]

    # make reports updated with overlap column
    merged_cnv = update_report(cnv, intersection, 'cnv')
    merged_sv = update_report(sv, intersection, 'sv')

    return collapse_overlaps(merged_cnv, 'cnv'), collapse_overlaps(merged_sv, 'sv')


def collapse_overlaps(report_df, variant_type):
    # if a variant (e.g. CNV) overlaps with muliple variants,
    # (e.g. multiple SVs), collapse and create one record
    # with a comma-separated list of overlapping variants
    overlap_col = 'SV_overlap' if variant_type == 'cnv' else 'CNV_overlap'
    start_col = 'START' if variant_type == 'cnv' else 'POS'
    # concatenate overlapping SVs
    merged_overlaps = report_df.groupby(['CHROM', start_col, 'END', 'SVTYPE'])[
        overlap_col].apply(','.join).reset_index()
    report_df = report_df.drop(overlap_col, axis=1)
    # merge report variants with merged overlaps df to replace overlap
    # annotations with list of overlaps, then drop duplicates
    report_df = report_df.merge(merged_overlaps, how='left', on=[
                                'CHROM', start_col, 'END', 'SVTYPE']).drop_duplicates()

    return report_df


def update_report(report_df, intersection, variant_type):
    # join overlap dataframe with report dataframe to annotate overlaps
    if variant_type == 'cnv':
        merged = report_df.merge(intersection, how='left', left_on=['CHROM', 'START', 'END', 'SVTYPE'], right_on=[
                                 'CNV_CHROM', 'CNV_START', 'CNV_END', 'CNV_SVTYPE'], validate='1:m')
        overlap = []
        overlap = ['{}:{}-{}'.format(row['SV_CHROM'], int(row['SV_START']), int(row['SV_END']))
                   if row['SV_START'] == row['SV_START'] else '.' for index, row in merged.iterrows()]
        merged['SV_overlap'] = overlap
    else:
        merged = report_df.merge(intersection, how='left', left_on=['CHROM', 'POS', 'END', 'SVTYPE'], right_on=[
                                 'SV_CHROM', 'SV_START', 'SV_END', 'SV_SVTYPE'], validate='1:m')
        overlap = []
        overlap = ['{}:{}-{}'.format(row['CNV_CHROM'], int(row['CNV_START']), int(row['CNV_END']))
                   if row['CNV_START'] == row['CNV_START'] else '.' for index, row in merged.iterrows()]
        merged['CNV_overlap'] = overlap
    merged = merged.drop(columns=['CNV_CHROM', 'CNV_START', 'CNV_END',
                                  'CNV_SVTYPE', 'SV_CHROM', 'SV_START', 'SV_END', 'SV_SVTYPE'])
    return merged


def reorder_columns(report_df):
    # Takes final annotated report and rearranges columns so overlap columns
    # is to the right of the 'SVTYPE' columns
    cols = report_df.columns.tolist()
    cols = cols[0:4] + [cols[-1]] + cols[4:-1]
    report_df = report_df[cols]
    return(report_df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Adds column to SV and CNV reports indicating if an SV is observed in the CNV report and vice versa')
    parser.add_argument('-sv', help='SV report (tsv file)', required=True)
    parser.add_argument('-cnv', help='CNV report (tsv file)', required=True)
    args = parser.parse_args()

    sv_file = args.sv.strip('tsv')
    cnv_file = args.cnv.strip('tsv')
    sv = pd.read_csv(args.sv, sep='\t')
    cnv = pd.read_csv(args.cnv, sep='\t')

    # determine overlaps between CNVs and SVs
    overlaps = intersect_variants(cnv, sv)

    # reorder columns
    cnv_w_overlap = reorder_columns(overlaps[0])
    sv_w_overlap = reorder_columns(overlaps[1])

    # export reports with overlap annotation
    cnv_w_overlap.to_csv('{}withSVoverlaps.tsv'.format(
        cnv_file), sep='\t', index=False, na_rep='nan')
    sv_w_overlap.to_csv('{}withCNVoverlaps.tsv'.format(
        sv_file), sep='\t', index=False, na_rep='nan')
