import pandas as pd, sys
import xlsxwriter

'''
format the outputs from ANNOVAR table_annovar.pl with OMIM and GNOMAD
merged.rare.EHdn.annovar.omim.hg19_multianno.txt
merged.rare.EHdn.annovar.gnomad.hg19_multianno.txt
'''

gnomad_file = sys.argv[1]
omim_file = sys.argv[2]
outfile = sys.argv[3]




omim = pd.read_csv(omim_file, sep="\t", header=[0,1])
omim.columns = [j if i.startswith("Otherinfo") else i for i, j in omim.columns ]
omim.rename(columns=lambda x: x.replace(".refGene",""), inplace=True)
omim_col = ['omim_inheritance', 'omim_phenotype' ]

gnomad = pd.read_csv(gnomad_file, sep="\t", header=[0,1])
gnomad.columns = [j if i.startswith("Otherinfo") else i for i, j in gnomad.columns ]
gnomad.rename(columns=lambda x: x.replace(".refGene",""), inplace=True)
gnomad_col = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func', 'Gene', 'motif', 'outliers', 'chr#start#end', 'size', 'a1000g_freq_perc', 'GeneDetail','ExonicFunc', 'AAChange', 'transcript', 'oe_lof', 'oe_lof_upper', 'oe_mis', 'oe_mis_upper', 'pLI', 'pRec']

annot_rare = pd.merge(gnomad[gnomad_col],omim[omim_col],right_index=True, left_index=True)
if not ".xlsx" in outfile:
    outfile = outfile + ".xlsx"
annot_rare.to_excel(outfile, index=False)
