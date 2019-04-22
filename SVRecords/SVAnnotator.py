import numpy as np
import pandas as pd
import sqlite3
import os
import subprocess
from pybedtools import BedTool
from collections import defaultdict
from enum import Enum, auto

class SVTYPE(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name
    
    DEL = auto()
    DUP = auto()
    INS = auto()
    INV = auto()
    IDP = auto()

class SVAnnotator:
    def __init__(self, exon_bed, hgmd_db, hpo, exac, omim, biomart):
        self.final_gene_ref_cols = ['BioMart Ensembl Gene ID', 'BioMart Associated Gene Name', 'HPO Features', \
       'OMIM Phenotypes', 'OMIM Inheritance', 'ExAC syn_z', 'ExAC mis_z', \
       'ExAC lof_z', 'ExAC pLI', 'HGMD disease', 'HGMD tag', 'HGMD descr', 'HGMD JOURNAL_DETAILS', 'HGMD SVTYPE']

        self.make_gene_ref_df(biomart)
        print('Annotating genes with HPO terms')
        self.annotate_hpo(hpo)
        print('Annotating genes with OMIM phenotypes and inheritance patterns')
        self.annotate_omim(omim)
        print('Annotating genes with ExAC transcript probabilities')
        self.annotate_exac(exac)
        self.annotate_hgmd(hgmd_db)

        # strip columns not used for annotation
        # technical note: drop duplicates before setting index - doing the reverse order will drop all duplicated columns instead of keeping one copy
        self.gene_ref_df = self.gene_ref_df[self.final_gene_ref_cols].drop_duplicates(keep='first').set_index('BioMart Ensembl Gene ID').astype(str)

    def set_column_values(self, df, annotation_dict, column_name):
        for interval, data in annotation_dict.items():
            df.loc[interval, column_name] = data

    def append_prefix_to_columns(self, df, prefix):
        #used to avoid collision of column names between different dataframes
        df.rename(columns={col: "%s %s" % (prefix, col) for col in df.columns}, inplace=True)

    def make_hgnc_dict(self, hgnc):
        hgnc_dict = {}
        with open(hgnc) as f:
            for line in f.readlines()[1:]:
                fields = line.strip('\n').split("\t")
                approved_symbol = fields[1].upper()
                prev_symbols = fields[4].split(', ')
                synonymous_symbols = fields[5].split(', ')

                for sym in prev_symbols:
                    hgnc_dict[sym.upper()] = approved_symbol
                
                for sym in synonymous_symbols:
                    hgnc_dict[sym.upper()] = approved_symbol
        
        self.hgnc_dict = hgnc_dict

    def gene2hgnc(self, col):
        def translate(cell):
            if isinstance(cell, str):
                if cell in self.hgnc_dict:
                    #print('individual sym match: %s %s' % (cell, self.hgnc_dict[cell]))
                    return self.hgnc_dict[cell]
                return cell
                #return self.hgnc_dict[cell] if cell in self.hgnc_dict else cell
            elif isinstance(cell, list):
                for sym in cell:
                    if sym in self.hgnc_dict: 
                        # if len(cell) > 1: print('list match: %s %s' % (sym, self.hgnc_dict[sym]))
                        # else: print('single element match: %s %s' % (sym, self.hgnc_dict[sym]))
                        return self.hgnc_dict[sym]
                return cell[0] # translation to hgnc gene has failed, just picked the first gene name then
            else:
                ValueError('DataFrame member is not of type str or list. %s' % str(cell))

        return col.apply(translate)

    def calc_exons_spanned(self, sample_df, exon_bed):
        exon_counts = defaultdict(int) #incrementing dict and setting value in df is faster than incrementing values in the df
        exon_ref = BedTool(exon_bed)
        sample_bedtool = BedTool(list(sample_df.reset_index()[['CHROM', 'POS', 'END', 'SVTYPE']].values))

        for interval in sample_bedtool.intersect(exon_ref, wa=True):
            exon_counts[(str(interval.chrom), str(interval.start), str(interval.stop),  str(interval[3]))] += 1
        
        count_df = pd.Series(exon_counts).to_frame().astype(str)
        count_df.index.names = ['CHROM', 'POS', 'END', 'SVTYPE']
        count_df.columns = ['EXONS_SPANNED']

        return sample_df.join(count_df).fillna(value={'EXONS_SPANNED': 0})

    def annotate_hgmd(self, hgmd):
        def get_hgmd_df():
            conn = sqlite3.connect(hgmd)

            gros_del = pd.read_sql_query('''
                SELECT del.DISEASE, del.TAG, del.DESCR, del.gene,
                printf('%s:%s:%s:%s', del.JOURNAL, del.AUTHOR, del.YEAR, del.PMID) AS JOURNAL_DETAILS,
                ALLGENES.HGNCID
                FROM GROSDEL as del 
                LEFT JOIN ALLGENES ON ALLGENES.GENE=del.GENE;
            ''', conn)

            gros_ins = pd.read_sql_query('''
                SELECT ins.DISEASE, ins.TAG, ins.DESCR, ins.gene, 
                printf('%s:%s:%s:%s', ins.JOURNAL, ins.AUTHOR, ins.YEAR, ins.PMID) AS JOURNAL_DETAILS,
                ALLGENES.HGNCID
                FROM GROSINS as ins 
                LEFT JOIN ALLGENES ON ALLGENES.GENE=ins.GENE
                WHERE ins.type='I';
            ''', conn)

            gros_dup = pd.read_sql_query('''
                SELECT ins.DISEASE, ins.TAG, ins.DESCR, ins.gene, 
                printf('%s:%s:%s:%s', ins.JOURNAL, ins.AUTHOR, ins.YEAR, ins.PMID) AS JOURNAL_DETAILS,
                ALLGENES.HGNCID
                FROM GROSINS as ins 
                LEFT JOIN ALLGENES ON ALLGENES.GENE=ins.GENE
                WHERE ins.type='D';
            ''', conn)
            
            conn.close()

            return gros_del, gros_ins, gros_dup

        def groupby_genes(df):
            df['hgncID'] = df['hgncID'].astype(str)
            #df['omimid'] = df['omimid'].astype(str)
            #df['gene'] = self.gene2hgnc(df['gene'])
            df = df.groupby(by='gene', as_index=False).agg(lambda x: "%s" % ', '.join(x))
            df['hgncID'] = df['hgncID'].apply(lambda col: col.split(', ')[0])
            #df['omimid'] = df['omimid'].apply(lambda col: col.split(', ')[0])
            return df
        
        matching_fields = {
            # "HGMD hgncID" : "BioMart HGNC ID(s)",
            "HGMD gene" : "BioMart Associated Gene Name",
        }

        hgmd_sv_df = pd.DataFrame()
        gros_del, gros_ins, gros_dup = get_hgmd_df()

        gros_del = groupby_genes(gros_del)
        gros_ins = groupby_genes(gros_ins)
        gros_dup = groupby_genes(gros_dup)

        for df in [gros_del, gros_ins, gros_dup]:
            df['gene'] = df['gene'].apply(lambda symbol: symbol.upper())
            self.append_prefix_to_columns(df, 'HGMD')

        gros_del['HGMD SVTYPE'] = SVTYPE.DEL.value
        gros_ins['HGMD SVTYPE'] = SVTYPE.INS.value
        gros_dup['HGMD SVTYPE'] = SVTYPE.DUP.value

        hgmd_sv_df = pd.concat([gros_del, gros_ins, gros_dup], ignore_index=True, sort=False).drop_duplicates().astype(str)

        self.gene_ref_df = self.prioritized_annotation(self.gene_ref_df, hgmd_sv_df, matching_fields)

    def prioritized_annotation(self, gene_ref_df, annotation_df, matched_fields):
        matched_rows = []

        ann_match_columns = list(matched_fields.keys())
        ref_match_columns = list(matched_fields.values())
        gene_ref_df_matching_cols = gene_ref_df[ref_match_columns].drop_duplicates()

        annotation_df = annotation_df.drop_duplicates()

        #join on equivalent fields
        for ann_field, ref_field in matched_fields.items():
            matched = gene_ref_df_matching_cols.join(annotation_df.set_index(ann_field), on=ref_field, how='inner')
            matched_rows.append(matched)

        #drop columns from annotation_df used for joining
        for table in matched_rows:
            for field in ann_match_columns:
                try:
                    table.drop(columns=[field, ], axis=0, inplace=True)
                    #print("Dropped %s" % field)
                except KeyError:
                    pass #tried to drop column which was joined on

        merged_df = pd.concat(matched_rows, ignore_index=True, sort=False).drop_duplicates().set_index(ref_match_columns).dropna(how='all')

        #add the remaining fields to the reference dataframe
        return self.gene_ref_df.join(merged_df, on=ref_match_columns, how='left').drop_duplicates()

    def annotate_hpo(self, hpo):
        matching_fields = {'HPO Gene ID': 'BioMart Ensembl Gene ID',
        # 'HPO Gene symbol': 'BioMart Associated Gene Name'
        }

        hpo_df = pd.read_csv(hpo, sep='\t')
        hpo_df.columns = hpo_df.columns.str.strip()
        # hpo_df = hpo_df[['Gene ID', 'Gene symbol', 'Features']]
        hpo_df = hpo_df[['Gene ID', 'Features']]
        hpo_df['Features'] = hpo_df['Features'].apply(lambda features: features.replace('; ', ', '))
        # hpo_df['Gene symbol'] = hpo_df['Gene symbol'].apply(lambda symbol: symbol.upper())

        hpo_df = hpo_df.astype(str)
        self.append_prefix_to_columns(hpo_df, "HPO")

        self.gene_ref_df = self.prioritized_annotation(self.gene_ref_df, hpo_df, matching_fields)
        #self.gene_ref_df.to_csv("hpo_ann.tsv", sep="\t")
    
    def annotate_omim(self, omim):
        omim_inheritance_codes = {"Autosomal dominant":"AD", \
        "Autosomal recessive":"AR", \
        "X-linked dominant":"XLD", \
        "X-linked recessive":"XLR", \
        "Y-linked dominant":"YLD", \
        "Y-linked recessive":"YLR", \
        "X-linked":"XL", \
        "Y-linked":"YL"}
        omim_inheritance_codes = {key.lower():value for key,value in omim_inheritance_codes.items()} # for case insensitive text search

        def process_OMIM_phenotype(phenotype):
            '''
                omim phenotype example:
                
                {Epilepsy, generalized, with febrile seizures plus, type 5, susceptibility to}, 613060 (3), Autosomal dominant; 
                {Epilepsy, idiopathic generalized, 10}, 613060 (3), Autosomal dominant; 
                {Epilepsy, juvenile myoclonic, susceptibility to}, 613060 (3), Autosomal dominant
            '''

            inheritance = []

            if pd.isnull(phenotype): return phenotype
            else:
                for p in phenotype.split('; '):
                    multiple_inheritance = [code for description, code in omim_inheritance_codes.items() if description.lower() in p.lower()]
                    if multiple_inheritance: inheritance.append('&'.join(multiple_inheritance))

                return ', '.join(inheritance)

        matching_fields = {'OMIM Ensembl Gene ID': 'BioMart Ensembl Gene ID'}

        # OMIM adds comments to their CSV file. These comments start with '#' character and are present in the header and footer of the file.
        omim_df = pd.read_csv(omim, sep='\t', header=3, skipfooter=61, engine='python')
        omim_df.columns = omim_df.columns.str.replace('#','')
        omim_df.columns = omim_df.columns.str.strip()

        omim_df = omim_df[['Mim Number', 'Ensembl Gene ID', 'Phenotypes']]
        omim_df = omim_df[pd.notnull(omim_df['Phenotypes'])] #drop all nan phenotype columns
        omim_df['Inheritance'] = omim_df['Phenotypes'].apply(lambda col: process_OMIM_phenotype(col))
        omim_df = omim_df.astype(str).groupby('Ensembl Gene ID', as_index=False).agg({'Phenotypes' : ' & '.join, 'Mim Number' : ' & '.join, 'Inheritance' : ' & '.join,})
        self.append_prefix_to_columns(omim_df, "OMIM")
        self.gene_ref_df = self.prioritized_annotation(self.gene_ref_df, omim_df, matching_fields)
        # self.gene_ref_df.to_csv("omim_ann.tsv", sep="\t")

    def annotate_exac(self, exac):
        matching_fields = {
            'ExAC gene' : 'BioMart Associated Gene Name',
        }

        exac_df = pd.read_csv(exac, sep='\t')
        exac_df.columns = exac_df.columns.str.strip()
        exac_df['transcript'] = exac_df['transcript'].apply(lambda transcript_id: transcript_id.split('.')[0])
        exac_df = exac_df[['gene', 'syn_z', 'mis_z', 'lof_z', 'pLI']]

        exac_df = exac_df.astype(str)
        self.append_prefix_to_columns(exac_df, "ExAC")

        self.gene_ref_df = self.prioritized_annotation(self.gene_ref_df, exac_df, matching_fields)
        #self.gene_ref_df.to_csv("exac_ann.tsv", sep="\t")

    def annotate_gnomad(self, gnomad, sv_record, reciprocal_overlap=0.9):
        gnomad_cols = ['CHROM',	'START', 'END', 'NAME',	'SVTYPE', 'AN',	'AC', 'AF', 'N_HOMREF',	'N_HET', 'N_HOMALT', 'FREQ_HOMREF',	'FREQ_HET',	'FREQ_HOMALT',	'POPMAX_AF']
        gnomad_ann_cols = ['gnomAD_SVTYPE', 'gnomAD_AN', 'gnomAD_AC', 'gnomAD_AF', 'gnomAD_N_HOMREF', 'gnomAD_N_HET', 'gnomAD_N_HOMALT', 'gnomAD_FREQ_HOMREF', 'gnomAD_FREQ_HET', 'gnomAD_FREQ_HOMALT', 'gnomAD_POPMAX_AF']
        gnomad_df = pd.read_csv(gnomad, sep='\t', dtype='str').astype(str)
        gnomad_df.columns = gnomad_df.columns.str.replace('#', '')
        gnomad_df.columns = gnomad_df.columns.str.strip()
        gnomad_df = gnomad_df[gnomad_cols]
        gnomad_bed = BedTool(gnomad_df.itertuples(index=False))

        sample_sv = sv_record.make_ref_bedtool()
        sv_record.df[gnomad_ann_cols] = pd.DataFrame([["NA" for field in gnomad_ann_cols]], index=sv_record.df.index)

        for ann in sample_sv.intersect(gnomad_bed, wa=True, wb=True, F=reciprocal_overlap, f=reciprocal_overlap):

            samp_chr, samp_start, samp_end, samp_svtype = ann[0:4]
            gnomad_chr, gnomad_start, gnomad_end, gnomad_id, gnomad_svtype = ann[4:9]

            if gnomad_svtype == samp_svtype:
                sv_record.df.loc[(samp_chr, samp_start, samp_end, samp_svtype), gnomad_ann_cols] = ann[8:20]
            # elif gnomad_svtype == 'MCNV':
            #     print("passing MCNV")

    def annotsv(self, sample_df):
        '''
            Handles DGV, DDD annotations
        '''
        all_sv_bed_name = "all_sv.bed"
        annotated = "./{}.annotated.tsv".format(all_sv_bed_name)
        sample_df.reset_index()[['CHROM', 'POS', 'END', 'SVTYPE']].to_csv(all_sv_bed_name, index=False, sep='\t')
        subprocess.call("$ANNOTSV/bin/AnnotSV -SVinputFile {} -SVinputInfo 1 -outputFile {}".format(all_sv_bed_name, annotated), shell=True)

        annotsv_df = pd.read_csv(annotated, sep='\t').astype(str)
        annotsv_df = annotsv_df.loc[annotsv_df['AnnotSV type'] == 'full']
        annotsv_df = annotsv_df.rename(columns={annotsv_df.columns[0]: 'CHROM', annotsv_df.columns[1]: 'POS', annotsv_df.columns[2]: 'END', annotsv_df.columns[3]:'SVTYPE'}).set_index(keys=['CHROM', 'POS', 'END', 'SVTYPE'])
        annotsv_df = annotsv_df[annotsv_df.columns[13:22]] # all dgv and ddd columns
        sample_df = sample_df.join(annotsv_df)

        os.remove(all_sv_bed_name)
        os.remove(annotated)

        return sample_df

    def make_gene_ref_df(self, biomart):
        df = pd.read_csv(biomart, sep='\t')
        # df = df[['Ensembl Gene ID', 'Ensembl Transcript ID', 'Associated Gene Name', 'HGNC ID(s)', 'MIM Gene Accession']].drop_duplicates()
        df = df[['Ensembl Gene ID', 'Associated Gene Name',]].drop_duplicates()
        df['Associated Gene Name'] = df['Associated Gene Name'].apply(lambda symbol: symbol.upper()) #make all gene symbols a single case to increase match rate with other dataframes
        df = df.astype(str)
        self.append_prefix_to_columns(df, "BioMart")
        
        self.gene_ref_df = df
    
    def annotate_genes(self, sample_df, gene_col):
        # extract genes from sample_df, create a new dataframe where each row only has a single ensemble id and interval info
        gene_df = sample_df.apply(lambda x: pd.Series(x[gene_col]),axis=1).stack().reset_index(level=4, drop=True)
        gene_df.name = gene_col
        gene_df = gene_df.to_frame().drop_duplicates().astype(str)

        # annotate passed in ensemble gene id's using the generated reference dataframe
        gene_df = gene_df.join(self.gene_ref_df, on=gene_col).drop_duplicates().reset_index()
        # parse out hgmd gene annotations who's SVTYPE does not match up with the sample's
        gene_df = gene_df[ (gene_df['SVTYPE'] == gene_df['HGMD SVTYPE']) | (pd.isnull(gene_df['HGMD SVTYPE'])) | (gene_df['HGMD SVTYPE'] == 'nan') ].set_index(['CHROM', 'POS', 'END', 'SVTYPE'])

        # aggregate all annotation columns within the same sv interval
        gene_df[gene_df.columns] = gene_df.groupby(gene_df.index).agg(list)[gene_df.columns]
        # parse out and replace nan values with "na" string
        gene_df[gene_df.columns] = gene_df[gene_df.columns].applymap(lambda cell: [str(item) for item in cell])
        gene_df[gene_df.columns] = gene_df[gene_df.columns].applymap(lambda cell: ["na"] if all("nan" == item.lower() for item in cell) else cell)
        gene_df[gene_df.columns] = gene_df[gene_df.columns].applymap(lambda cell: ["na" if "nan" == item.lower() else item for item in cell])
        gene_df[gene_df.columns] = gene_df[gene_df.columns].applymap(lambda cell: ', '.join(cell))
        gene_df = gene_df[gene_df.columns].drop_duplicates()

        # annotate the passed in dataframe
        sample_df = sample_df.drop(gene_col, axis=1).join(gene_df)
        return sample_df
