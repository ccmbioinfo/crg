import argparse
import numpy as np
import pandas as pd
from pybedtools import BedTool
import re


def get_sample_names(sv_report, variant_type):
    colnames = sv_report.columns
    samples = [col for col in colnames if 'SV_DETAILS' in col]
    if variant_type == 'cnv' and ('WGS' in samples[0] or 'RGS' in samples[0]):
        #sometimes tcag cnv reports have '_WGS' or '_RGS' in sample name
        samples = [sample.split('_')[0] + '_' + sample.split('_')[1] + '_' + sample.split('_')[2] for sample in samples]
    else:    
        samples = [sample.split('_')[0] + '_' + sample.split('_')[1] for sample in samples]
        

    return samples


def get_sample_sv_detail_names(sample_id, variant_type):
    sample_col = '{}_SV_DETAILS'.format(sample_id)
    details = 'SV_DETAILS' if variant_type == 'sv' else 'CNV_DETAILS'
    return sample_col, details


def variant_details_to_df(variants, sample_id, variant_type):
    # takes <sample>_SV_DETAILS column and turns it into a dataframe
    # that can be converted into a BedTool for intersection
    sample_col, details = get_sample_sv_detail_names(sample_id, variant_type)
    variants = variants[[sample_col]].copy()
    variants = variants[variants[sample_col] != '.']
    # for unfiltered reports, there can be multiple comma-separated
    # SVs in the SV_DETAILS column; take the first one
    for index, row in variants[sample_col].items():
        if len(row.split(':')) != 3:
            variants[sample_col][index] = row.split(",")[0]
    # split columns
    variants[['CHR', 'COORD', 'SVTYPE']] = variants[sample_col].str.split(':', expand=True)
    variants[['START', 'END']] = variants['COORD'].str.split('-', expand=True)
    variants = variants.drop('COORD', axis=1)
    variants = variants[['CHR', 'START', 'END', 'SVTYPE', sample_col]]
    return variants


def intersect_variants(cnv, sv):
    # intersects CNVs and SVs and returns dataframes of intersections

    # convert report dataframes to bedtools
    cnv_bed = BedTool.from_dataframe(cnv)
    sv_bed = BedTool.from_dataframe(sv)

    bed_cols = ['CNV_CHROM', 'CNV_START', 'CNV_END', 'CNV_SVTYPE', 'CNV_DETAILS',
                'SV_CHROM', 'SV_START', 'SV_END', 'SV_SVTYPE', 'SV_DETAILS']
    # intersect SVs and CNVs with 50% reciprocal overlap
    intersection_rec50 = cnv_bed.intersect(sv_bed, wa=True, wb=True, F=0.5,
                                           f=0.5).saveas('intersection_rec50.bed')
    intersection_rec50 = pd.read_csv('intersection_rec50.bed', sep='\t', names=bed_cols)
    # make sure CNV and SV are same variant type
    intersection_rec50 = intersection_rec50[intersection_rec50['CNV_SVTYPE']
                                            == intersection_rec50['SV_SVTYPE']]
    intersection_rec50['INTERSECTION'] = ['rec50']*len(intersection_rec50)                                        
    # intersect SVs and CNVs with any overlap
    intersection_any = cnv_bed.intersect(sv_bed, wa=True, wb=True).saveas('intersection_any.bed')
    intersection_any = pd.read_csv('intersection_any.bed', sep='\t', names=bed_cols)
    # make sure CNV and SV are same variant type
    intersection_any = intersection_any[intersection_any['CNV_SVTYPE']
                                            == intersection_any['SV_SVTYPE']]    
    intersection_any['INTERSECTION'] = ['any']*len(intersection_any)  

    intersection_all = pd.concat([intersection_rec50, intersection_any], ignore_index=True)     
    intersection_all = intersection_all.groupby(bed_cols)['INTERSECTION'].apply(','.join).reset_index()                                                               

    return intersection_all


def update_report(report_df, intersection, variant_type, sample_id):
    # join overlap dataframe with report dataframe to annotate overlaps

    sample_col, details = get_sample_sv_detail_names(sample_id, variant_type)
    merged = report_df.merge(intersection, how='left', left_on=sample_col, right_on=details)
    start_col = 'START' if variant_type == 'cnv' else 'POS'
    merged = merged.drop_duplicates(subset=['CHROM', start_col, 'END', 'SVTYPE'])

    return merged


def add_tag(variant_details, intersection, variant_type):
    # if overlap between CNV and SV exists, add tag to sample _SV_DETAILS
    # CNV_50/SV_50 if 50% reciprocal overlap, CNV_any/SV_any if any intersections

    # if neither field is empty, add tag

    variant_overlap = 'CNV' if variant_type == 'sv' else 'SV'

    if (variant_details == variant_details) and (intersection == intersection):
        intersection = '50' if 'rec50' in intersection else 'any'
        variant_details = variant_details + ":{}_{}".format(variant_overlap, intersection)
    else:
        pass
    return variant_details


def apply_add_tag(updated_report_df, variant_type, sample_id):
    # adds tags to report dataframe
    sample_col, details = get_sample_sv_detail_names(sample_id, variant_type)
    updated_report_df[sample_col] = updated_report_df.apply(
        lambda row: add_tag(
            variant_details=row[sample_col],
            intersection=row['INTERSECTION'],
            variant_type=variant_type), axis=1)
    return(updated_report_df)


def trim_columns(tagged_report_df):
    # removes columns from report df that were added as a result of the merge
    # between the report and the bed file of CNV-SV intersections
    tagged_report_df = tagged_report_df.drop(columns=['CNV_CHROM', 'CNV_START', 'CNV_END', 'CNV_SVTYPE', 'CNV_DETAILS',
                                                      'SV_CHROM', 'SV_START', 'SV_END', 'SV_SVTYPE', 'SV_DETAILS', 'INTERSECTION'])
    return tagged_report_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Adds column to SV and CNV reports indicating if an SV is observed in the CNV report and vice versa')
    parser.add_argument('-sv', help='SV report (tsv file)', required=True)
    parser.add_argument('-cnv', help='CNV report (tsv file)', required=True)
    args = parser.parse_args()

    # get file names
    sv_file = args.sv.strip('tsv')
    cnv_file = args.cnv.strip('tsv')

    # read reports into dataframes
    sv = pd.read_csv(args.sv, sep='\t')
    cnv = pd.read_csv(args.cnv, sep='\t')

    # get sample names
    sv_sample_names = get_sample_names(sv, 'sv')
    cnv_sample_names = get_sample_names(cnv, 'cnv')

    for sv_sample, cnv_sample in zip(sv_sample_names, cnv_sample_names):
        # convert sample SV_DETAILS into dataframe
        cnvs = variant_details_to_df(cnv, cnv_sample, 'cnv')
        svs = variant_details_to_df(sv, sv_sample, 'sv')

        # get intersection between SVs and CNVs
        intersection = intersect_variants(cnvs, svs)
        intersection.to_csv('{}_intersection.tsv'.format(sv_sample), sep='\t')

        # tag sample SV_DETAILS if overlap is present
        cnv = update_report(cnv, intersection, 'cnv', cnv_sample)
        cnv = apply_add_tag(cnv, 'cnv', cnv_sample)
        cnv = trim_columns(cnv)

        sv = update_report(sv, intersection, 'sv', sv_sample)
        sv = apply_add_tag(sv, 'sv', sv_sample)
        sv = trim_columns(sv)

    if 'unfiltered' in sv_file.split('.'):
        sv.to_csv('{}withCNVoverlaps.tsv'.format(sv_file), index=False, na_rep='nan', sep='\t')
    else:
        cnv.to_csv('{}withSVoverlaps.tsv'.format(cnv_file), index=False, na_rep='nan', sep='\t')
        sv.to_csv('{}withCNVoverlaps.tsv'.format(sv_file), index=False, na_rep='nan', sep='\t')

