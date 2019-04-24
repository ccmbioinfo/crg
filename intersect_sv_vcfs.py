import argparse
from SVRecords import SVGrouper, SVAnnotator

def main(exon_bed, hgmd_db, hpo, exac, omim, biomart, gnomad, outfile_name, vcfs):
    SVScore_fields = ['variants/SVSCORESUM', 'variants/SVSCOREMAX', 'variants/SVSCORETOP5', 'variants/SVSCORETOP10', 'variants/SVSCOREMEAN',]

    print("Grouping like structural variants ...")
    sv_records = SVGrouper(vcfs, ann_fields=SVScore_fields)

    print('Annotating structural variants ...')
    ann_records = SVAnnotator(exon_bed, hgmd_db, hpo, exac, omim, biomart)
    sv_records.df = ann_records.annotate_genes(sv_records.df, "Ensembl Gene ID")
    # sv_records.df = ann_records.annotsv(sv_records.df)
    # sv_records.df = ann_records.calc_exons_spanned(sv_records.df, exon_bed)
    # ann_records.annotate_gnomad(gnomad, sv_records)
    ann_records.add_decipher_link(sv_records.df)

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
    parser.add_argument('-o', help='Output file name e.g. -o 180.sv.family.tsv', required=True, type=str)
    args = parser.parse_args()

    if len(args.i) == 0:
        ValueError('Please enter the path to some vcf\'s following the -i flag')
    else:
        main(args.exon_bed, args.hgmd, args.hpo, args.exac, args.omim, args.biomart, args.gnomad, args.o, args.i)