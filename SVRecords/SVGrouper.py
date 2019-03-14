import allel
import numpy as np
import pandas as pd
import os
from pybedtools import BedTool

class SVGrouper:
    def __init__(self, vcfs):
        def list2string(col):
            return col.apply(lambda x: ', '.join(x) if isinstance(x, list) else x)

        all_sv, sample_list = self._merge_sv(vcfs)
        assert len(sample_list) == len(set(sample_list)), "Duplicate sample names among input vcf's detected: %s" % sample_list

        columns = ['CHROM', 'POS', 'END', 'SVTYPE', 'Ensembl Gene ID', 'N_SAMPLES']
        columns.extend(sample_list)
        columns.extend(["%s_SV_DETAILS" % s for s in sample_list])
        columns.extend(["%s_GENOTYPE" % s for s in sample_list])

        self.sample_list = sample_list
        self.df = pd.DataFrame(columns=columns)
        self.df.set_index(keys=['CHROM', 'POS', 'END', 'SVTYPE'], inplace=True)
        self._group_sv(all_sv)
        #list set typecast is to ensure that we get a list of unique ensemble identifiers
        self.df['Ensembl Gene ID'] = self.df['Ensembl Gene ID'].apply(lambda gene_list: list(set(gene_list.split(','))) if ',' in gene_list else gene_list)
        self.bedtool = self._make_ref_bedtool()

        for name in self.sample_list:
            self.df["%s_SV_DETAILS" % name] = list2string(self.df["%s_SV_DETAILS" % name])
            self.df["%s_GENOTYPE" % name] = list2string(self.df["%s_GENOTYPE" % name])

    def _merge_sv(self, vcf_paths):
        '''
            Merge all SV interval data from multiple vcf's in to a single BedTool instance

            Implementation:
                Use Panda's dataframe for some easy preprocessing, then create a BedTool from a tuple containing each row
        '''
        def split_Ensembl_ids(id_list):
            new_list = []
            for id in id_list:
                if '-' in id:
                    new_list.extend(id.split('-'))
                elif '&' in id:
                    new_list.extend(id.split('&'))
                else:
                    new_list.append(id)
            return new_list

        intervals = []
        sample_names = []

        for vcf_path in vcf_paths:
            vcf_dict = allel.read_vcf(vcf_path, ['samples', 'calldata/GT', 'variants/CHROM', 'variants/END', 'variants/POS', 'variants/SVTYPE', 'variants/ANN', ], numbers={'ANN': 1000}, transformers=allel.ANNTransformer()) #use read_vcf because genotype field is not picked up with vcf_to_dataframe
            
            # vcf_dict = allel.read_vcf(vcf_path, ['variants/CHROM', 'variants/END', 'variants/POS', 'variants/SVTYPE', 'samples', 'calldata/GT', 'variants/ANN_Gene_ID',\
            #     'variants/SVLEN', 'variants/ANN', 'variants/SVSCORESUM', 'variants/SVSCOREMAX', 'variants/SVSCORETOP5', 'variants/SVSCORETOP10', 'variants/SVSCOREMEAN', ], \
            #     numbers={'ANN': 1000}, transformers=allel.ANNTransformer()) #use read_vcf because genotype field is not picked up with vcf_to_dataframe

            assert len(vcf_dict['samples']) == 1, "%s contains 0 or more than 1 sample: %s" % (vcf_path, str(vcf_dict['samples']))
            name = vcf_dict.pop('samples')[0]
            sample_names.append(name)

            #discard all other ANN fields
            for key in list(vcf_dict.keys()):
                if "ANN" in key and (key != "variants/ANN_Gene_ID"):
                    vcf_dict.pop(key)

            # remove empty strings, split on delimited characters, then join using comma
            vcf_dict['variants/ANN_Gene_ID'] = [list(filter(None, ann)) for ann in vcf_dict['variants/ANN_Gene_ID']] #by default, specifying numbers=1000 creates 1000 elements, with most being empty
            vcf_dict['variants/ANN_Gene_ID'] = [split_Ensembl_ids(id_list) if any('&' in id for id in id_list) or any('-' in id for id in id_list) else id_list for id_list in vcf_dict['variants/ANN_Gene_ID']]
            vcf_dict['variants/ANN_Gene_ID'] = [','.join(id_list) if isinstance(id_list, list) else id_list for id_list in vcf_dict['variants/ANN_Gene_ID']]

            vcf_dict['calldata/GT'] = np.array(['HET' if 0 in gt and 1 in gt else 'HOM' for gt in vcf_dict.pop('calldata/GT')])

            df = pd.DataFrame(vcf_dict)
            df['samples'] = name
            df = df[['variants/CHROM', 'variants/POS', 'variants/END', 'variants/SVTYPE', 'calldata/GT', 'variants/ANN_Gene_ID', 'samples']] #order of first 3 columns matters for creation of BedTool

            intervals.extend(df.itertuples(index=False))

        return BedTool(intervals), sample_names

    def _group_sv(self, bedtool):
        already_grouped_intervals = set()

        for l in bedtool.intersect(bedtool, wa=True, wb=True, F=0.5, f=0.5):

            ref_chr, ref_start, ref_end, ref_svtype, ref_gt, ref_genes, ref_name, \
            samp_chr, samp_start, samp_end, samp_svtype, samp_gt, samp_genes, samp_name = l

            ref_interval = (ref_chr, ref_start, ref_end, ref_svtype)
            samp_interval = (samp_chr, samp_start, samp_end, samp_svtype, samp_gt, samp_name)

            if (samp_interval not in already_grouped_intervals) and (ref_svtype == samp_svtype):
                self._add_interval(ref_interval, ref_genes, samp_interval)
                already_grouped_intervals.add(samp_interval)
        
        self.df.sort_index(inplace=True)

    def _add_interval(self, ref_interval, ref_genes, samp_interval):
        samp_chr, samp_start, samp_end, samp_svtype, samp_gt, samp_name = samp_interval

        #Get reference to row
        if ref_interval not in self.df.index:
            #make new row
            self.df.loc[ref_interval, :] = np.nan
            row = self.df.loc[ref_interval, :]
            row['N_SAMPLES'] = 0
            row['Ensembl Gene ID'] = ref_genes

            for name in self.sample_list:
                row[name] = 0
                row["%s_SV_DETAILS" % name] = []
                row["%s_GENOTYPE" % name] = []
        else:
            row = self.df.loc[ref_interval, :]

        #Set values for row
        if row[samp_name] == 0:
            row['N_SAMPLES'] += 1
            row[samp_name] = 1

        row["%s_SV_DETAILS" % samp_name].append('{}:{}-{}:{}:{}'.format(samp_chr, samp_start, samp_end, samp_svtype, samp_gt))
        row["%s_GENOTYPE" % samp_name].append(samp_gt)
    
    def write(self, outfile_name):
        self.df.to_csv(outfile_name, sep='\t', encoding='utf-8')

    def _make_ref_bedtool(self):
        return BedTool(list(self.df.index.values))