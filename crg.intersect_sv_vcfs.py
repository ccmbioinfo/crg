import argparse
from SVRecords import SVGrouper, SVAnnotator

def main(exon_bed, hgmd_db, hpo, exac, omim, biomart, gnomad, sv_counts, outfile_name, vcfs):
    SVScore_cols = ['variants/SVLEN', 'variants/SVSCORESUM', 'variants/SVSCOREMAX', 'variants/SVSCORETOP5', 'variants/SVSCORETOP10', 'variants/SVSCOREMEAN',]
    HPO_cols = [ "N_UNIQUE_HPO_TERMS", "HPO Features", "N_GENES_IN_HPO", "Genes in HPO" ]

    print("Grouping like structural variants ...")
    sv_records = SVGrouper(vcfs, ann_fields=SVScore_cols)
    sample_cols = [ col for col in sv_records.df.columns if col != "Ensembl Gene ID" ]
    sample_genotype_cols = [col for col in sample_cols if col.endswith('_GENOTYPE')]
    # sv_records.df.to_csv('grouped.tsv', sep='\t')

    print('Annotating structural variants ...')
    ann_records = SVAnnotator(exon_bed, hgmd_db, hpo, exac, omim, biomart)
    sv_records.df = ann_records.annotate_genes(sv_records.df, "Ensembl Gene ID")

    for sv_count in sv_counts:
        sv_records.df = ann_records.annotate_counts(sv_count, sv_records, prefix=sv_count)
    
    sv_records.df = ann_records.annotsv(sv_records.df)
    sv_records.df = ann_records.calc_exons_spanned(sv_records.df, exon_bed)
    sv_records.df = ann_records.annotate_gnomad(gnomad, sv_records)
    sv_records.df = ann_records.annotate_hgmd(hgmd_db, sv_records.df)
    ann_records.add_decipher_link(sv_records.df)

    if not set(HPO_cols).issubset(set(sv_records.df.columns)):
        for col in HPO_cols:
            sv_records.df[col] = "na"

    # format and rearrange the columns
    # sv_records.df = sv_records.df[ [col for col in sample_cols if col not in SVScore_cols and col not in sample_genotype_cols] + \
    # [ "EXONS_SPANNED", "Ensembl Gene ID", "BioMart Associated Gene Name", ] + \
    # [ "N_UNIQUE_HPO_TERMS", "HPO Features", "N_GENES_IN_HPO", "Genes in HPO" ] + \
    # [ "N_GENES_IN_OMIM", "Genes in OMIM", "OMIM Phenotypes", "OMIM Mim Number", "OMIM Inheritance" ] + \
    # [ "DGV_LOSS_IDs", "DGV_LOSS_n_samples_with_SV", "DGV_LOSS_n_samples_tested", "DGV_LOSS_Frequency", "DDD_SV", "DDD_DUP_n_samples_with_SV", "DDD_DUP_Frequency", "DDD_DEL_n_samples_with_SV", "DDD_DEL_Frequency" ] + \
    # [ "DECIPHER_LINK" ] + \
    # [ "ExAC syn_z", "ExAC mis_z", "ExAC lof_z", "ExAC pLI" ] + \
    # [ "Genes in HGMD", "HGMD disease", "HGMD tag", "HGMD descr", "HGMD JOURNAL_DETAILS" ] + \
    # [ "gnomAD_SV", "gnomAD_AN", "gnomAD_AC", "gnomAD_AF", "gnomAD_N_HOMREF", "gnomAD_N_HET", "gnomAD_N_HOMALT", "gnomAD_FREQ_HOMREF", "gnomAD_FREQ_HET", "gnomAD_FREQ_HOMALT", "gnomAD_POPMAX_AF"] + \
    # sample_genotype_cols + \
    # SVScore_cols ] 
    sv_records.df.columns = sv_records.df.columns.str.replace('variants/','')

    print('Writing results to file ...')
    sv_records.write(outfile_name)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generates a structural variant report in CSV format for clincal research')
    parser.add_argument('-i', type=str, nargs='+', help='VCF files containing structural variants, e.g. -i 180_230.vcf 180_231.vcf', required=True)
    parser.add_argument('-exon_bed', help='BED file containing fixed exon positions', required=True)
    parser.add_argument('-hgmd', help='HGMD Pro database file', required=True, type=str)
    parser.add_argument('-hpo', help='Tab delimited file containing gene names and HPO terms', type=str)
    parser.add_argument('-exac', help='ExAC tab delimited file containing gene names and scores', type=str, required=True)
    parser.add_argument('-omim', help='OMIM tab delimited file containing gene names and scores', type=str, required=True)
    parser.add_argument('-biomart', help='TSV file from BiomaRt containing Ensemble gene ID, transcript ID, gene name, MIM gene id, HGNC id, EntrezGene ID', type=str, required=True)
    parser.add_argument('-gnomad', help='BED file from Gnomad containing structural variant coordinates and frequencies across populations', type=str, required=True)
    parser.add_argument('-sv_counts', nargs='+', help='List of BED files containing structural variants and their frequencies. Can be used to annotate with various populations and variant callers', required=False)
    parser.add_argument('-overlap', help='Recipricol overlap to group a structural variant by', type=float, default=0.5)
    parser.add_argument('-o', help='Output file name e.g. -o 180.sv.family.tsv', required=True, type=str)
    args = parser.parse_args()

    if len(args.i) == 0:
        ValueError('Please enter the path to some vcf\'s following the -i flag')
    else:
        main(args.exon_bed, args.hgmd, args.hpo, args.exac, args.omim, args.biomart, args.gnomad, args.sv_counts, args.o, args.i)