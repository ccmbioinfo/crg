import os
import subprocess
import sqlite3
from pybedtools import BedTool

class GeneAnnotations:
	def __init__(self, gene_name):
		self.gene_name = gene_name

		#DDD fields, not yet implemented
		self.DDD_status = ""
		self.DDD_mode = ""
		self.DDD_consequence = ""
		self.DDD_disease = ""
		self.DDD_pmids = ""

		self.is_in_hgmd = False
		#HGMD GROSS mutation fields
		self.hgmd_disease = []
		self.hgmd_tag = []
		self.hgmd_description = []
		self.hgmd_comments = []
		self.hgmd_journal = []
		self.hgmd_author = []
		self.hgmd_year = []
		self.hgmd_pmid = []

		#OMIM fields
		self.mim_num = ""
		self.mim_inheritance = ""
		self.mim_phenotype = ""

		#pLI scores
		self.synz = ""
		self.misz = ""
		self.pli = ""

	def add_hgmd_anno(self, hgmd_annots):
		for row in hgmd_annots:
			disease, tag, description, comments, journal, author, year, pmid = row
			self.hgmd_disease.append(disease)
			self.hgmd_tag.append(tag)
			self.hgmd_description.append(description)
			self.hgmd_comments.append(comments)
			self.hgmd_journal.append(journal)
			self.hgmd_author.append(author)
			self.hgmd_year.append(year)
			self.hgmd_pmid.append(pmid)

	def group_journal_fields(self):
		journals = [':'.join([j, a, y, p]) for j, a, y, p in zip(self.hgmd_journal, self.hgmd_author, self.hgmd_year, self.hgmd_pmid)]
		return ";".join(journals) if any(journals) else ""
			
class StructuralVariant:
	def __init__(self, chr, start, end, svtype, genotype, svlen, svsum, svmax, svtop5, svtop10, svmean, dgv):
		self.chr = chr
		self.start = start
		self.end = end
		self.genotype = genotype

		#SVSCORES fields
		self.svlen = svlen
		self.svtype = svtype
		self.svmax = svmax
		self.svsum = svsum
		self.svtop5 = svtop5
		self.svtop10 = svtop10
		self.svmean = svmean

		#Annotated by TCAG
		self.dgv = dgv

		# GeneAnnotation values 
		self.genes = {}

		# Custom implementation
		self.exons_spanned = 0

		#AnnotSV - DGV Fields
		self.dgv_gain_id = "NA"
		self.dgv_gain_n_samples_with_sv = "NA"
		self.dgv_gain_n_samples_tested = "NA"
		self.dgv_gain_freq = "NA"

		self.dgv_loss_id = "NA"
		self.dgv_loss_n_samples_with_sv = "NA"
		self.dgv_loss_n_samples_tested = "NA"
		self.dgv_loss_freq = "NA"

		#AnnotSV - DDD fields
		self.ddd_sv = ""
		self.ddd_dup_n_samples_with_sv = ""
		self.ddd_dup_freq = ""
		self.ddd_del_n_samples_with_sv = ""
		self.ddd_del_freq = ""

	def make_decipher_link(self):
		return '=hyperlink("https://decipher.sanger.ac.uk/browser#q/{}:{}-{}")'.format(self.chr, self.start, self.end)

	def key(self):
		return (self.chr, self.start, self.end)

	def make_interval_string(self):
		return '{}:{}-{}:{}:{}'.format(self.chr, self.start, self.end, self.svtype, self.genotype)
	
	def add_gene(self, gene_name):
		if gene_name not in self.genes.keys():
			self.genes[gene_name] = GeneAnnotations(gene_name)
			return self.genes[gene_name]
		else:
			raise ValueError('Attempted to add %s twice to %s' % (gene_name, self.make_interval_string()))

	def make_gene_list(self):
		return ';'.join([gene for gene in self.genes])

	def make_column_from_list(self, column_data):
		def isListEmpty(inList):
			#Returns true if inList contains some combination of empty nested lists and/or empty strings
			#Taken and modified from: https://stackoverflow.com/questions/1593564/python-how-to-check-if-a-nested-list-is-essentially-empty
			if isinstance(inList, list): # is a list
				return all( map(isListEmpty, inList) )
			elif isinstance(inList, str) and not inList: # is an empty string
				return True
			elif isinstance(inList, str) and inList: # is a non-empty string
				return False
			else:
				raise ValueError("Heterogenous list containing data elements other than list or str detected! Raw data: %s" % str(inList))

		if isListEmpty(column_data):
			return "NA"
		elif all(isinstance(i, list) for i in column_data): #nested list
			return ';'.join(['|'.join([data if data else "NA" for data in gene_data]) if any(gene_data) \
			else "NA" \
			for gene_data in column_data])
		elif all(isinstance(i, str) for i in column_data): #all strs
			return ';'.join([data if data \
			else "NA" \
			for data in column_data])

	def make_hgmd_gene_list(self):
		return self.make_column_from_list([gene.gene_name for gene in self.genes.values() if gene.is_in_hgmd])

	def make_gene_mim_columns(self):
		#return mim_num, mim_inheritance, mim_phenotype
		return self.make_column_from_list([gene.mim_num for gene in self.genes.values()]), \
		self.make_column_from_list([gene.mim_inheritance for gene in self.genes.values()]), \
		self.make_column_from_list([gene.mim_phenotype for gene in self.genes.values()])

	def make_gene_hgmd_columns(self):
		#return disease, tag, description, comments, journal_info
		return self.make_column_from_list([gene.hgmd_disease for gene in self.genes.values()]), \
		self.make_column_from_list([gene.hgmd_tag for gene in self.genes.values()]), \
		self.make_column_from_list([gene.hgmd_description for gene in self.genes.values()]), \
		self.make_column_from_list([gene.hgmd_comments for gene in self.genes.values()]), \
		self.make_column_from_list([gene.group_journal_fields() for gene in self.genes.values()])

	def make_pli_columns(self):
		#return pli, misz, synz
		return self.make_column_from_list([ gene.pli for gene in self.genes.values() ]), \
		self.make_column_from_list([ gene.misz for gene in self.genes.values() ]), \
		self.make_column_from_list([ gene.synz for gene in self.genes.values() ])

	def make_ddd_columns(self):
		#return self.ddd_sv, self.ddd_dup_n_samples_with_sv, self.ddd_dup_freq, self.ddd_del_n_samples_with_sv, self.ddd_del_freq
		return [column if column else "NA" for column in [self.ddd_sv, self.ddd_dup_n_samples_with_sv, self.ddd_dup_freq, self.ddd_del_n_samples_with_sv, self.ddd_del_freq]]

class StructuralVariantRecords:
	'''
		Holds groupings of StructuralVariant.

		Groupings in the dict grouped_sv follows this structure:
			grouped_sv[ref_interval][samp_name] = [StructuralVariant]

		Annotated information on the ref_interval used to group sample intervals is stored in:
			all_ref_interval_data[(chr, start, end)] = StructuralVariant
	'''
	def __init__(self, _sample_list):
		self.grouped_sv = {}
		self.all_ref_interval_data = {}
		self.sample_list = _sample_list

	def add_interval(self, ref_interval, new_interval, column_data, samp_name):
		'''
			Adds new_interval to the dict under the key ref_interval
			Column data relating to a ref_interval is stored as well
		'''
		if ref_interval not in self.grouped_sv.keys():
			self.grouped_sv[ref_interval] = {}
			self.all_ref_interval_data[ref_interval] = column_data[ref_interval]

		if samp_name not in self.grouped_sv[ref_interval]:
			self.grouped_sv[ref_interval][samp_name] = []

		self.grouped_sv[ref_interval][samp_name].append(column_data[new_interval])

	def make_header(self):
		fields = ["CHR", "START", "END", "N_SAMPLES", "LIST", "N_GENES", "GENES", "LONGEST_SVTYPE", "SVLEN", "SVSCORE_MAX", \
		"SVSCORE_SUM", "SVSCORE_TOP5", "SVSCORE_TOP10", "SVSCORE_MEAN", "DGV", "EXONS_SPANNED", "DECIPHER_LINK", \
		"DGV_GAIN_IDs", "DGV_GAIN_n_samples_with_SV", "DGV_GAIN_n_samples_tested", "DGV_GAIN_Frequency", "DGV_LOSS_IDs", \
		"DGV_LOSS_n_samples_with_SV", "DGV_LOSS_n_samples_tested", "DGV_LOSS_Frequency", "DDD_SV", "DDD_DUP_n_samples_with_SV", \
		"DDD_DUP_Frequency", "DDD_DEL_n_samples_with_SV", "DDD_DEL_Frequency", "N_GENES_IN_OMIM", "OMIM_MIM_NUMBER", "OMIM_INHERITANCE", "OMIM_PHENOTYPE", \
		"synZ", "misZ", "pLI", "GENES_IN_HGMD", "HGMD_DISEASE", "HGMD_TAG", "HGMD_DESCRIPTION", "HGMD_COMMENTS", "HGMD_JOURNAL_INFO"]
		fields.extend(self.sample_list)
		fields.extend(["%s_SV_details" % s for s in self.sample_list])
		fields.extend(["%s_genotype" % s for s in self.sample_list])
		return "%s\n" % ",".join(fields)

	def make_bed(self, bed_name):
		'''
			Creates a bed file containing all reference intervals
		'''
		with open(bed_name, "w") as f:
			for sv in self.grouped_sv:
				f.write('{}\t{}\t{}\t{}\n'.format(sv[0], sv[1], sv[2], bed_name))

	def all_ref_BedTool(self):
		return BedTool([sv for sv in self.all_ref_interval_data.keys()])

	def write_results(self, outfile_name):
		'''
			A lot of string manipulation to generate the final CSV file line by line
		'''
		def make_sample_list_index(sample_list, interval):
			'''
				Makes strings for LIST and SAMPLENAME column
				e.g. 1;2;3;4;5;6;7 and 1,1,1,1,1,1,1
			'''
			return ";".join([str(i+1) for i, sample in enumerate(self.sample_list) if sample in interval]), \
			",".join(["1" if sample in interval else "0" for sample in sample_list])

		def make_sample_sv_details(sample_list, interval):
			'''
				Makes list for SAMPLENAME_details column
				e.g. 1:10334731-10334817:DEL;1:10334769-10334833:DUP
			'''
			return ','.join(["NA" if sample not in interval \
			else ";".join([variant.make_interval_string() for variant in interval[sample]]) \
			for sample in sample_list ])

		def make_sample_genotype_details(sample_list, interval):
			return ','.join(["NA" if sample not in interval \
			else "HOM" if "HOM" in [variant.genotype for variant in interval[sample]] \
			else "HET" \
			for sample in sample_list])

		def get_longest_svtype(interval):
			'''
				Returns SVTYPE of largest SV in a grouping
			'''
			def svlen(type_len):
				return abs(int(type_len[1]))

			svtype_and_svlen = [(variant.svtype, variant.svlen) for sample in interval.values() for variant in sample]
			svtype_and_svlen.sort(key=svlen)
			return svtype_and_svlen[-1][0] #largest SV is last value in sorted list, hence [-1]

		with open(outfile_name, "w") as out:
			out.write(self.make_header())

			for key in sorted(self.grouped_sv.keys()):

				chr, start, end = key
				ref = self.all_ref_interval_data[key]

				sample_list_index, sample_list_isthere_index = make_sample_list_index(self.sample_list, self.grouped_sv[key])
				n_samples = str(len(self.grouped_sv[key]))
				svtype = get_longest_svtype(self.grouped_sv[key])
				samp_sv_details = make_sample_sv_details(self.sample_list, self.grouped_sv[key])
				samp_genotype_details = make_sample_genotype_details(self.sample_list, self.grouped_sv[key])
				n_genes = str(len(ref.genes))
				n_mim_genes = str(len([ gene for gene in ref.genes.values() if gene.mim_num ]))
				ddd_sv, ddd_dup_n_samples_with_sv, ddd_dup_freq, ddd_del_n_samples_with_sv, ddd_del_freq = ref.make_ddd_columns()

				#GeneAnnotations
				synz, misz, pli = ref.make_pli_columns()
				hgmd_disease, hgmd_tag, hgmd_description, hgmd_comments, hgmd_journal_info = ref.make_gene_hgmd_columns()
				mim_num, mim_inheritance, mim_phenotype = ref.make_gene_mim_columns()
				genes = ref.make_gene_list()

				out_line = u'{},{},{},{}\n'.format(",".join([str(chr), str(start), str(end), n_samples, sample_list_index, n_genes, genes, svtype, \
				ref.svlen, ref.svmax, ref.svsum, ref.svtop5, ref.svtop10, ref.svmean, ref.dgv, str(ref.exons_spanned), ref.make_decipher_link(), \
				ref.dgv_gain_id, ref.dgv_gain_n_samples_with_sv, ref.dgv_gain_n_samples_tested, ref.dgv_gain_freq, ref.dgv_loss_id, \
				ref.dgv_loss_n_samples_with_sv, ref.dgv_loss_n_samples_tested, ref.dgv_loss_freq, ddd_sv, ddd_dup_n_samples_with_sv, \
				ddd_dup_freq, ddd_del_n_samples_with_sv, ddd_del_freq, n_mim_genes, mim_num, mim_inheritance, mim_phenotype, synz, misz, pli, \
				ref.make_hgmd_gene_list(), hgmd_disease, hgmd_tag, hgmd_description, hgmd_comments, hgmd_journal_info]), \
				sample_list_isthere_index, samp_sv_details, samp_genotype_details)

				out.write(out_line)

	def annotate(self, exon_bed, hgmd_db):
		def calc_exons_spanned(exon_bed):
			'''
				Populates the field: exons_spanned for all reference intervals
			'''
			exon_ref = BedTool(exon_bed)
			all_ref_sv = self.all_ref_BedTool()
			for interval in all_ref_sv.intersect(exon_ref, wa=True):
				self.all_ref_interval_data[str(interval.chrom), str(interval.start), str(interval.stop)].exons_spanned += 1

		def annotsv():
			'''
				Handles DGV, DDD and OMIM annotations
			'''
			all_sv_bed_name = "all_sv.bed"
			annotated = "./{}.annotated.tsv".format(all_sv_bed_name)
			self.make_bed(all_sv_bed_name)

			subprocess.call("$ANNOTSV/bin/AnnotSV -SVinputFile {} -SVinputInfo 1 -outputFile {}".format(all_sv_bed_name, annotated), shell=True)

			with open(annotated, "r") as f:
				next(f) #skip header

				for fields in f:
					field = fields.rstrip('\n').replace(',', ';').split('\t')
					field = [ "" if not f else str(f) for f in field ]
					sv = self.all_ref_interval_data[(field[0], field[1], field[2])]
					if field[4] == "full":
						#DGV Annotations
						sv.dgv_gain_id = field[13]
						sv.dgv_gain_n_samples_with_sv = field[14]
						sv.dgv_gain_n_samples_tested = field[15]
						sv.dgv_gain_freq = field[16]

						sv.dgv_loss_id = field[17]
						sv.dgv_loss_n_samples_with_sv = field[18]
						sv.dgv_loss_n_samples_tested = field[19]
						sv.dgv_loss_freq = field[20]

						#DDD Annotations
						sv.ddd_sv = field[21]
						sv.ddd_dup_n_samples_with_sv = field[22]
						sv.ddd_dup_freq = field[23]
						sv.ddd_del_n_samples_with_sv =field[24]
						sv.ddd_del_freq = field[25]
					elif field[4] == "split":
						gene_name = field[5]
						gene_anno = sv.genes[gene_name] if gene_name in sv.genes.keys() else sv.add_gene(gene_name)

						#OMIM Annotations
						gene_anno.mim_num, gene_anno.mim_phenotype, gene_anno.mim_inheritance = field[38], field[39], field[40]
						#pLI Annotations
						gene_anno.synz, gene_anno.misz, gene_anno.pli = field[34], field[35], field[36]

			os.remove(all_sv_bed_name)
			os.remove(annotated)

		def hgmd(db_path):
			def decode_rows(rows):
				return [ [ "" if field is None or not field else \
				str(field) if isinstance(field,int) else \
				field.replace(',', ';') \
				for field in row ] \
				for row in rows ]

			conn = sqlite3.connect(db_path)
			cur = conn.cursor()

			for sv in self.all_ref_interval_data.values():
				svtype = sv.svtype

				for gene_anno in sv.genes.values():
					gene_name = gene_anno.gene_name
					cur.execute('SELECT GENE FROM ALLGENES WHERE GENE=?', (gene_name, ))

					if cur.fetchone() is not None:
						gene_anno.is_in_hgmd = True
						# Look for solved cases in this gene involving structural variants 
						if svtype == "DEL":
							cur.execute('SELECT DISEASE, TAG, DESCR, COMMENTS, JOURNAL, AUTHOR, YEAR, PMID FROM GROSDEL WHERE GENE=?', (gene_name, ) )
							gene_anno.add_hgmd_anno(decode_rows(cur.fetchall()))
						elif svtype == "INS":
							cur.execute('SELECT DISEASE, TAG, DESCR, COMMENTS, JOURNAL, AUTHOR, YEAR, PMID FROM GROSINS WHERE GENE=? AND TYPE=?', (gene_name, 'I'))
							gene_anno.add_hgmd_anno(decode_rows(cur.fetchall()))
						elif svtype == "DUP":
							cur.execute('SELECT DISEASE, TAG, DESCR, COMMENTS, JOURNAL, AUTHOR, YEAR, PMID FROM GROSINS WHERE GENE=? AND TYPE=?', (gene_name, 'D'))
							gene_anno.add_hgmd_anno(decode_rows(cur.fetchall()))
					else:
						gene_anno.is_in_hgmd = False

			conn.close()

		print('Querrying HGMD ...')
		hgmd(hgmd_db)
		print('Running AnnotSV for DDD, DGV, OMIM, PLI columns')
		annotsv()
		print('Calculating exons spanned ...')
		calc_exons_spanned(exon_bed)