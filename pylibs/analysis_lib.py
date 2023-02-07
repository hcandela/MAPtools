from logging import exception
from pylibs.constants import *
import re
import sys
import os
from math import log
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import hypergeom
from scipy.spatial import distance
import pandas as pd
import gffpandas.gffpandas as gffpd
from Bio import SeqIO
from Bio.Seq import Seq
#from docopt import docopt

def cmHH(command_line):
	os.system(command_line)

def test_args(__doc__,arg):
	if arg['--input'] == None and arg['pipe'] == True:
		print(__doc__, end='')
		sys.exit()
	if arg['--input'] != None:
		inp_f = arg['--input']
		try:
			f = open(inp_f, 'r')
		except FileNotFoundError:
			print('Error: The input file {} does not exist'.format(inp_f))
			sys.exit()
	
	
	wd = os.getcwd()
	try:
		os.makedirs(wd+'/'+arg['--outdir'])
	except FileExistsError:
		print('#Warning: the output directory already exists')
		pass

	arg['--outdir'] = wd+'/'+arg['--outdir']+'/'
	if arg['--output'] != None:
		arg['--fileformat'] = '.'+arg['--fileformat']
		arg['--output'] = arg['--output'] + arg['--fileformat']
		arg['--output'] = check_save_an(arg, arg['--output'])
	else:
		arg['--output'] = None
	
	if arg['--fileformat'] == '.csv':
		arg['spacer'] = ','
	else:
		arg['spacer'] = '\t'
	arg['lim'] = 1e-90
	if 'mbs' in arg.keys():
		arg = check_mbs_args(arg)
	if 'qtl' in arg.keys():
		arg = check_qtl_args(arg)
	return arg

def check_save_an(arg, file_name):
	typ = arg['--fileformat']
	if os.path.isfile(arg['--outdir']+file_name):
		expand = 1
		while True:
			expand += 1
			nw_file_name = file_name.split(typ)[0] +'_'+ str(expand) + typ
			if os.path.isfile(arg['--outdir']+nw_file_name):
				continue
			else:
				file_name = nw_file_name
				return file_name
	else:
		return file_name

def check_mbs_args(arg):
	data = arg['--data'].split(',')
	if not 'R' in data:
		print('Error: You should include the recessive pool (--data R,X,X')
		sys.exit()
	if 'Pr' in data and 'Pd' in data:
		print('Error: You should include only one parental pool in data (--data D,R,Px or --data R,Px')
		sys.exit()
	arg['--min-depth'] = int(arg['--min-depth'])
	arg['--max-depth'] = int(arg['--max-depth'])
	arg['--min-ratio'] = int(arg['--min-ratio'])/100
	arg['--max-ratio'] = int(arg['--max-ratio'])/100
	if arg['--min-depth'] <= 0:
		arg['--min-depth'] = 1
	if arg['--max-depth'] <= arg['--min-depth']:
		print('Error: You should choose a correct interval of depths(--min_depth 20 --max-depth 120)')
		sys.exit()
	if arg['--min-ratio'] <= 0:
		arg['--min-ratio'] = 0
	if arg['--max-ratio'] <= arg['--min-ratio']:
		print('Error: You should choose a correct interval of frequencies(--min_ratio 15 --max-ratio 85)')
		sys.exit()
	return arg
def check_qtl_args(arg):
	data = arg['--data'].split(',')
	if len(data) < 2:
		print('Error: You should include a minimum of 2 pools in QTL-seq experiment (--data D,R,Px or --data D,R)')
		sys.exit()
	if 'Pr' in data and 'Pd' in data:
		print('Error: You should include only one parental pool in data (--data D,R,Px)')
		sys.exit()

	arg['--min-depth'] = int(arg['--min-depth'])
	arg['--max-depth'] = int(arg['--max-depth'])
	if arg['--min-depth'] <= 0:
		arg['--min-depth'] = 1
	if arg['--max-depth'] <= arg['--min-depth']:
		print('Error: You should choose a correct interval of depths(--min_depth 20 --max-depth 120)')
		sys.exit()
	return arg
def LF(entrada):
	result = 0
	if entrada == 0:
		return 0
	else:
		for i in range(entrada, 1, -1):
			if i >= 1:
				result += log(i)
			else:
				break
	return result

def LogFisher(a,b,c,d):
	'''
	Calculates and returns fisher statistic logarithm as float
	Arguments: a, b, c, d. these numbers come from the allele count
	'''
	logfisher = (LF(a+b) + LF(c+d) + LF(a+c) + LF(b+d) - LF(a) - LF(b) - LF(c) - LF(d) - LF(a+b+c+d)) / log(10)
	return logfisher

def Gstatic(a,b,c,d):
	row1 = a + b
	row2 = c + d
	col1 = a + c
	col2 = b + d
	tot = a+b+c+d
	n1 = (row1*col1)/tot
	n2 = (row1*col2)/tot
	n3 = (row2*col1)/tot
	n4 = (row2*col2)/tot
	counts = [a,b,c,d]
	nx = [n1,n2,n3,n4]
	g = 0
	for i in range(len(counts)):
		if counts[i] > 0:
			g += 2*(counts[i]*log(counts[i]/nx[i]))
	return g
	
def pvalor(a,b,c,d):
	table = np.array([[a, b], [c, d]])
	oddsr, pvalue = fisher_exact(table, alternative='two-sided')
	return pvalue

def mbs_calc(inp, arg):
	'''
	Takes the allele count and filter thresholds for:
		max and min coverage (max_dp, min_dp)
		max and min ratios to consider a position heterozigous (max_ratio, min_ratio)
	Performs the necessary calcs for representation
	Returns a list of the calcs with:
		[
		ratio of the recesive allele of the recesive pool,
		ratio of the recesive allele of the dominant pool,
		ratio of the most frequent allele of the recesive pool,
		LogFisher of the allele count
		Function that generates a peak in the maximum of the representation
		p-value
		log10 of p-value
		]
	'''
	#TODO - Change boost to calculate always 
	data = arg['--data'].split(',')
	inf_s = set(data)
	min_dp = arg['--min-depth']
	max_dp = arg['--max-depth']
	min_ratio = arg['--min-ratio']
	max_ratio = arg['--max-ratio']
	if len(inf_s) == 1 and arg['--no-ref'] == True:
		a, b = inp[0], inp[1]
		if a != 0 or b != 0:
			ratio3 = max(a,b)/(a+b)
			boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
			return [ratio3, boost]
	elif (len(inf_s) == 1 and arg['--no-ref'] == False) or (len(inf_s) == 2 and 'D' not in inf_s):
		a, b = inp[0], inp[1]
		if a != 0 or b != 0:
			ratio1 = b/(a+b)
			boost = 1/(arg['lim'] + abs(1 - 1/max(ratio1, 1-ratio1)))
			return [ratio1, boost]
	elif 'D' in inf_s and 'R' in inf_s:
		a, b, c, d = inp[0], inp[1], inp[2], inp[3]
		dom = a + b #allele count of the dominant pool
		rec = c + d #allele count of the recesive pool
		if dom <= max_dp and dom >= min_dp and (a/dom) <= max_ratio and (a/dom) >= min_ratio: #filters for positions with enought coverage and heterozigous on the dominant pool
			ratio3 = max(c,d)/(c+d)
			resultado = LogFisher(a,b,c,d)
			pva = pvalor(a,b,c,d)
			pva10 = (log(pva)/log(10))
			if arg['--no-ref'] == True:
				boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
				return [ratio3,resultado,boost,pva,pva10]
			else:
				ratio1 = b/(a+b)
				ratio2 = d/(c+d)
				boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
				return [ratio1,ratio2,ratio3,resultado,boost,pva,pva10]

def qtl_calc(inp, arg):
	data = arg['--data'].split(',')
	inf_s = set(data)
	no_ref = arg['--no-ref']
	min_dp = arg['--min-depth']
	max_dp = arg['--max-depth']
	a,b,c,d = inp[0], inp[1], inp[2], inp[3]
	dom = a + b
	rec = c + d
	dp_total = dom + rec
	if dp_total <= max_dp and dp_total >= min_dp and rec > 0 and dom > 0:
		v1 = (c/rec, d/rec)
		v2 = (a/dom, b/dom)
		ed = distance.euclidean(v1,v2)
		g = Gstatic(a,b,c,d)
		pva = pvalor(a,b,c,d)
		pva10 = (log(pva)/log(10))
		if no_ref == True:
			return [ed,g,pva,pva10]
		else:
			ratio1 = b/(a+b)
			ratio2 = d/(c+d)
			delta = ratio1 - ratio2
			return [ratio1,ratio2,delta,ed,g,pva,pva10]

def new_line(fsal, arg, first, fields, al_count, calcs):
	data = arg['--data'].split(',')
	inf_s = set(data)
	res = fields + al_count + calcs
	spacer = arg['spacer']
	if first:
		if 'mbs' in arg.keys():
			if len(inf_s) == 1 or 'D' not in inf_s:
				if arg['--no-ref'] == True:
					header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','MAX_SNPidx2','BOOST']
				else:
					header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','SNPidx1','BOOST']

			elif 'D' in inf_s and 'R' in inf_s and arg['--no-ref'] == False:
				header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','SNPidx1','SNPidx2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE']

			elif 'D' in inf_s and 'R' in inf_s and arg['--no-ref'] == True:
				header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE']

		elif 'qtl' in arg.keys():
			if arg['--no-ref'] == True:
				header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','ED','G','PVALUE','log10PVALUE']

			if arg['--no-ref'] == False:
				header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','SNPidx1','SNPidx2','DELTA','ED','G','PVALUE','log10PVALUE']

		first = False
		n_head = spacer.join(header)+'\n'
		write_line(n_head, fsal)

	n_line = spacer.join(str(field) for field in res) + '\n'
	write_line(n_line, fsal)
	return first


def write_line(n_line, fsal):
	if fsal == False:
		print(n_line, end='')
	else:
		fsal.write(n_line)

def vcf_line_parser2(line, arg):
	'''
	This script reads the output from the bcftools call line by line and processes it
	The line is splitted into list z[]. Values of the matrix are:
		z[0] = chromosome
		z[1] = positionarguments
		z[2] = ID
		z[3] = Reference allele
		z[4] = Variant allele
		z[5] = Quality
		z[6] = Filter
		z[7] = Info line. Read Bcftools documentation for more information
		z[8] = Format. Contains the description of the format for the nexts fields
		z[9] to z[x] = Genotype and allele count for bam 1, 2, 3...x
	'''
	z = line.split('\t')
	alt = z[4] #alternative allele
	if alt == '.':
		#if the alternative allele is '.' indicates that there is no variation in the aligned sequences, so that position is discarded
		return 0,0
	else:
		form = z[8].split(':')
		GT_index = form.index('GT')
		AD_index = form.index('AD')
		data = arg['--data'].split(',')
		res = dict()
		inf = dict()
		c = 9
		for i in range(len(data)):		#dict inf contain the name of the pool and the column where is it in input 
			inf[data[i]] = c			#{'R': 9, 'D':10, 'Pr':11}
			c += 1
		inf_s = set(data)
		if len(inf_s) == 1:
			data1 = z[9].split(':')			#split data from R pool
			pool = dict()
			AD = data1[AD_index].split(',') #obtain AD data from R pool
			ref = z[3]
			alt = z[4].split(',')			#alt could be more than one nt
			pool[ref] = int(AD[0])			#generete a dict for the pool {'nt':AD1, 'nt':AD2 }
			for i in range(len(alt)):
				pool[alt[i]] = int(AD[i+1])
			res[data[0]] = pool				#Save the pool dict naming the pool {'R':{'nt1':AD1, 'nt2':AD2 }}
		elif len(inf_s) >= 2:
			if inf_s == {'D','R'} or inf_s == {'D','R','Pr'} or inf_s == {'D','R','Pd'}:
				data1 = z[inf['D']].split(':')	#save info from each pool
				data2 = z[inf['R']].split(':')
				GT1 = data1[GT_index].split('/')
				GT2 = data2[GT_index].split('/')
				gt = set(GT1+GT2)				#compare the genotypes
				if len(gt) == 1 or '.' in gt: #if all the genotypes are equal the line is discarded
					return 0,0
				for p,c in inf.items():		#for each pool and column
					ref = z[3]
					alt = z[4].split(',')
					AD = z[c].split(':')[AD_index].split(',')
					pool = {ref:int(AD[0])}
					for i in range(len(alt)):
						pool[alt[i]] = int(AD[i+1])
					res[p] = pool	#{'R':{'nt1':AD1, 'nt2':AD2 }, 'D':{'nt1':AD1, 'nt2':AD2}...}
			elif inf_s == {'R','Pr'} or inf_s == {'R','Pd'}:
				ref = z[3]
				alt = z[4].split(',')
				for p,c in inf.items():
					AD = z[c].split(':')[AD_index].split(',')
					pool = {ref:int(AD[0])}
					for i in range(len(alt)):
						pool[alt[i]] = int(AD[i+1])
					res[p] = pool
		return [z[0], z[1], z[3]],res #returns chromosome, position, reference allele, and the data for each bam)

def normalize(pools, REF, arg, r_min=0.03):
	data = arg['--data'].split(',')
	inf_s = set(data)
	ref = arg['--ref']
	n_ref = arg['--no-ref']
	rec = pools['R']
	alle = [a for a,v in rec.items()]
	REF_idx = alle.index(REF)
	if len(alle) > 2:
		return
	if len(inf_s) == 1:
		if ref == 'D' or n_ref == True:
			a = rec[alle[REF_idx]]
			alle.pop(REF_idx)
			b = rec[alle[0]]
			return [REF, alle[0], a, b]
		elif ref == 'R':
			b = rec[alle[REF_idx]]
			alle.pop(REF_idx)
			a = rec[alle[0]]
			return [alle[0],REF,a, b]

	elif len(inf_s) == 2:
		if n_ref == True:
			dom = pools['D']
			dom_list = [count for al,count in dom.items()]
			rec_list = [count for al,count in rec.items()]
			return[alle[0],alle[1],dom_list[0],dom_list[1],rec_list[0],rec_list[1]]
		elif inf_s == {'R','D'} and ref == 'D':
			dom = pools['D']
			a = dom[alle[REF_idx]]
			c = rec[alle[REF_idx]]
			alle.pop(REF_idx)
			b = dom[alle[0]]
			d = rec[alle[0]]
			return [REF, alle[0], a, b, c, d]
		elif inf_s == {'R','D'} and ref == 'R':
			dom = pools['D']
			b = dom[alle[REF_idx]]
			d = rec[alle[REF_idx]]
			alle.pop(REF_idx)
			a = dom[alle[0]]
			c = rec[alle[0]]
			return [alle[0], REF, a, b, c, d]
		elif 'Pd' in inf_s:
			p_dom = pools['Pd']
			p_al = sorted(p_dom, key=lambda key: p_dom[key], reverse=True) #ordered from higher to lower
			if p_dom[p_al[0]] > 0 and p_dom[p_al[1]]/(p_dom[p_al[0]] + p_dom[p_al[1]]) < r_min or len(p_al) > 2:
				a = rec[p_al[0]]
				b = rec[p_al[1]]
				return [p_al[0], p_al[1], a, b]
		elif 'Pr' in inf_s:
			p_rec = pools['Pr']
			p_al = sorted(p_rec, key=lambda key: p_rec[key], reverse=True)
			if p_rec[p_al[0]] > 0 and p_rec[p_al[1]]/(p_rec[p_al[0]] + p_rec[p_al[1]]) < r_min or len(p_al) > 2:
				a = rec[p_al[0]]
				b = rec[p_al[1]]
				return [p_al[1], p_al[0], a, b]
	elif len(inf_s) == 3:
		if 'Pr' in inf_s:
			parental = pools['Pr']
		elif 'Pd' in inf_s:
			parental = pools['Pd']
		p_al = sorted(parental, key=lambda key: parental[key], reverse=True)
		dom = pools['D']
		if parental[p_al[0]] > 0 and parental[p_al[1]]/(parental[p_al[0]] + parental[p_al[1]]) < r_min or len(p_al) > 2:
			if 'Pd' in inf_s:
				a = dom[p_al[0]]
				b = dom[p_al[1]]
				c = rec[p_al[0]]
				d = rec[p_al[1]]
				return[p_al[0], p_al[1], a, b, c, d]
			elif 'Pr' in inf_s:
				a = dom[p_al[1]]
				b = dom[p_al[0]]
				c = rec[p_al[1]]
				d = rec[p_al[0]]
				return[p_al[1], p_al[0], a, b, c, d]

def test_arg_ann(__doc__,arg):
	if arg['--input'] == None and arg['pipe'] == True:
		print(__doc__,end='')
		sys.exit()
	if arg['--input'] != None:
		inp_f = arg['--input']
		try:
			f = open(inp_f, 'r')
		except FileNotFoundError:
			print('Error: The input file {} does not exit'.format(inp_f))
			sys.exit()
	wd = os.getcwd()
	try:
		os.makedirs(wd+'/'+arg['--outdir'])
	except FileExistsError:
		print('#Warning: the output directory already exists')
		pass

	
	arg['--outdir'] = wd+'/'+arg['--outdir']+'/'
	if arg['--output'] != None:
		arg['--fileformat'] = '.txt'
		arg['--output'] = arg['--output'] + arg['--fileformat']
		arg['--output'] = check_save_an(arg, arg['--output'])
	else:
		arg['--output'] = None
	try:
		gff = open(wd+'/'+arg['--gff'])
		gff.close()
		arg['--gff'] = wd+'/'+arg['--gff']
	except FileNotFoundError:
		print('The gff file does not exist.')
		sys.exit()

	try:
		ref = open(wd+'/'+arg['--reference'])
		ref.close()
		arg['--reference'] = wd+'/'+arg['--reference']
	except FileNotFoundError:
		print('The reference file does not exist.')
		sys.exit()
	
	Tchrom,Treg = arg['--region'].split(':')
	Treg = Treg.split('-')
	arg['--region'] = [Tchrom,int(Treg[0]),int(Treg[1])]
	data = arg['--data'].split(',')
	arg['spacer'] = '\t'
	arg['mbs'] = True
	arg['lim'] = 1e-90
	#chroms = gff['seq_id'].unique()
	#gff_chrom = chroms[int(Tchrom) - 1]
	#arg['gff_chrom'] = gff_chrom
	
	if arg['--mutagen'] == None:
		arg['--mutagen'] = 'EMS'
	if not 'R' in data:
		print('Error: You should include the recessive pool (--data')
		sys.exit()
	return arg

def filter_region(line, arg):
	z = line.split('\t')
	chrom = z[0]
	pos = int(z[1])
	if chrom == arg['--region'][0] and pos >= arg['--region'][1] and pos <= arg['--region'][2]:
		return line
	else:
		return
def load_gff(arg):
	gff_path = arg['--gff']
	gff = {'chromosome':[],'gene': [], 'mRNA':[], 'five_prime_UTR':[], 'exon':[], 'CDS':[], 'three_prime_UTR':[],\
		 'miRNA':[], 'tRNA':[], 'ncRNA':[], 'ncRNA_gene':[], 'lnc_RNA':[], 'snoRNA':[], 'snRNA':[],\
		 'rRNA':[], 'pseudogene':[], 'pseudogenic_transcript':[], 'pre_miRNA':[], 'SRP_RNA':[]}
	with open(gff_path, 'r') as handle:
		for line in handle:
			if line.startswith('#'):
				continue
			else:
				line = line.rstrip()
				data = line.split('\t')
				seqid = data[0]
				if seqid == arg['--region'][0]: #TODO Complete condition to cath only structures in our region
					type_ = data[2]
					start = data[3]
					end = data[4]
					strand = data[6]
					phase = data[7]
					attributes = data[8].split(';')
					dict_att = {'ID':'.','Parent':'.','Name':'.'}
					for att in attributes:
						dict_att[att.split('=')[0]] = att.split('=')[1]
					#dict_att = {dict_att[att.split('=')[0]] :att.split('=')[1] for att in attributes}
					if type_ in gff.keys():

						gff[type_].append([seqid, type_, int(start), int(end), strand, phase, dict_att['ID'], dict_att['Parent'], dict_att['Name']])
		gff['gene'] = gff['gene'] + gff['ncRNA_gene'] + gff['pseudogene']
		gff['gene'].sort(key=lambda row:row[2])
		gff['mRNA'].sort(key=lambda row:(row[6], row[2]))
		gff['exon'].sort(key=lambda row:(row[7], row[2]))
		gff['CDS'].sort(key=lambda row:(row[7], row[2]))
		gff['five_prime_UTR'].sort(key=lambda row:(row[7], row[2]))
		gff['three_prime_UTR'].sort(key=lambda row:(row[7], row[2]))
		ordered = ['gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']
		for type_ in gff:
			if type_ not in ordered:
				gff[type_].sort(key=lambda row:row[2])
	arg['gff'] = gff

				
def  find_row(rows, pos, start, end):
	right = start
	current = end
	while right < current:
		m = (right + current)//2
		beg = rows[m][2]
		#end = rows[m][3]
		if pos < beg: current = m
		else:
			right = m + 1
	
	left = start
	current = end
	while left < current:
		m = (left + current)//2
		end = rows[m][3]
		if pos > end: left = m + 1
		else:
			current = m
	return left, right-1

def  find_row_name(rows, pos, key):
	right = 0
	current = len(rows)
	while right < current:
		m = (right + current)//2
		beg = rows[m][key]
		#end = rows[m][3]
		if pos < beg: current = m
		else:
			right = m + 1
	
	left = 0
	current = len(rows)
	while left < current:
		m = (left + current)//2
		end = rows[m][key]
		if pos > end: left = m + 1
		else:
			current = m
	return left, right-1

def load_reference(df,arg):
	chrom = df['#CHROM'].unique()[0]
	flag = True
	for sequ in SeqIO.parse(arg['--reference'], 'fasta'):
		if sequ.id == chrom:
			arg['ref'] = sequ
			flag = False
			break
	if flag == True:
		print('Chromosome not found in FASTA file.')
		sys.exit()


def filter_mut(arg, al):
	REF = al[0]
	ALT = al[1]
	if arg['--mutagen'] == 'EMS':
		if (REF == 'G' and ALT == 'A') or (REF == 'C' and ALT == 'T'):
			return True
		else:
			return False

def ann_calc(inp, arg):
	'''
	Takes the allele count and filter thresholds for:
		max and min coverage (max_dp, min_dp)
		max and min ratios to consider a position heterozigous (max_ratio, min_ratio)
	Performs the necessary calcs for representation
	Returns a list of the calcs with:
		[
		LogFisher of the allele count
		Function that generates a peak in the maximum of the representation
		p-value
		log10 of p-value
		]
	'''
	data = arg['--data'].split(',')
	inf_s = set(data)
	if len(inf_s) == 1 and arg['--no-ref'] == True:
		a, b = inp[0], inp[1]
		if a != 0 or b != 0:
			ratio3 = max(a,b)/(a+b)
			boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
			return [ratio3, boost]
	elif (len(inf_s) == 1 and arg['--no-ref'] == False) or (len(inf_s) == 2 and 'D' not in inf_s):
		a, b = inp[0], inp[1]
		if a != 0 or b != 0:
			ratio1 = b/(a+b)
			boost = 1/(arg['lim'] + abs(1 - 1/max(ratio1, 1-ratio1)))
			return [ratio1, boost]
	elif 'D' in inf_s and 'R' in inf_s:
		a, b, c, d = inp[0], inp[1], inp[2], inp[3]
		dom = a + b #allele count of the dominant pool
		rec = c + d #allele count of the recesive pool
		ratio3 = max(c,d)/(c+d)
		resultado = LogFisher(a,b,c,d)
		pva = pvalor(a,b,c,d)
		pva10 = (log(pva)/log(10))
		if arg['--no-ref'] == True:
			boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
			return [ratio3,resultado,boost,pva,pva10]
		else:
			ratio1 = b/(a+b)
			ratio2 = d/(c+d)
			boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
			return [ratio1,ratio2,ratio3,resultado,boost,pva,pva10]

def create_df(arg):
	data = arg['--data'].split(',')
	inf_s = set(data)
	if len(inf_s) == 1 or 'D' not in inf_s:
		if arg['--no-ref'] == True:
			header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','MAX_SNPidx2','BOOST']
		else:
			header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','SNPidx1','BOOST']
	elif 'D' in inf_s and 'R' in inf_s and arg['--no-ref'] == False:
		header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','SNPidx1','SNPidx2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE']
	elif 'D' in inf_s and 'R' in inf_s and arg['--no-ref'] == True:
		header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE']
	arg['header'] = header
	df = pd.DataFrame(columns=header)
	return df

def new_df_line(df, arg, fields, al_count, calcs):
	data = arg['--data'].split(',')
	inf_s = set(data)
	res = fields + al_count + calcs
	n_line = {arg['header'][i]:res[i] for i in range(len(res))}
	df_nline = pd.DataFrame.from_records([n_line])
	df = pd.concat([df, df_nline])
	return df


def check_mutation2(row, arg):
	gff = arg['gff']
	pos = int(row['POS'])
	b,e = find_row(gff['gene'], pos, 0, len(gff['gene']))
	if gff['gene'][b:e+1]:
		print(pos,'is genic')
		for gene in gff['gene'][b:e+1]:
			gene_id = gene[6]
			b1,e1 = find_row_name(gff['mRNA'], gene_id, 7)
			b2,e2 = find_row(gff['mRNA'], pos, b1, e1+1)
			for mRNA in gff['mRNA'][b2:e2+1]:
				print(mRNA)
				mRNA_id = mRNA[6]
				b3,e3 = find_row_name(gff['exon'], mRNA_id, 7)
				b4,e4 = find_row(gff['exon'], pos, b3, e3+1)
				if gff['exon'][b4:e4+1]:
					#for exon in gff['exon'][b4:e4+1]:
					#	print(exon)
					b5,e5 = find_row_name(gff['CDS'], mRNA_id, 7)
					b6,e6 = find_row(gff['CDS'], pos, b5, e5+1)
					if gff['CDS'][b6:e6+1]:
						for cds in gff['CDS'][b6:e6+1]:
							print(cds)
							beg, end = codon_coords(pos,cds[2],cds[3],cds[4],int(cds[5])) #target pos, cds start, cds end, cds strand, cds phase
							if beg < cds[2]:
								print("  Codon is truncated at exon/intron boundary at: ",cds[2])
								print("  Possible splice site mutation. Check adjacent exon for a possible aa substitution")
							elif end > cds[3]:
								print("  Codon is truncated at exon/intron boundary at: ",cds[3])
								print("  Possible splice site mutation. Check adjacent exon for a possible aa substitution")
							else:
								print(' Codon fully contained exon')
								codon = arg['ref'].seq[beg-1:end]
								before = codon[:pos-beg]+row['REF']+codon[pos-beg+1:]
								after = codon[:pos-beg]+row['ALT']+codon[pos-beg+1:]
								print(before, after)
								find_effect(before,after,cds[4])

					b7,e7 = find_row_name(gff['five_prime_UTR'], mRNA_id, 7)
					b8,e8 = find_row(gff['five_prime_UTR'], pos, b7, e7+1)
					b9,e9 = find_row_name(gff['three_prime_UTR'], mRNA_id, 7)
					b10,e10 = find_row(gff['three_prime_UTR'], pos, b9, e9+1)
					if gff['five_prime_UTR'][b8:e8+1]:
						five =  gff['five_prime_UTR'][b8:e8+1][0]
						print('is 5\' UTR')
						print(five)
						if five[4] == '+':
							distfive1 = five[3] - pos
							print(' Dist:', distfive1)
						elif five[4] == '-':
							distfive1 = pos - five[2]
							print(' Dist:', distfive1)
					elif gff['three_prime_UTR'][b10:e10+1]:
						three =  gff['three_prime_UTR'][b10:e10+1][0]
						print('is 3\' UTR')
						print(three)
		
				else:
					print('is intron')
					print(gff['exon'][b4])
					print(gff['exon'][e4])
					dis1 = gff['exon'][b4][2] - pos
					dis2 = pos - gff['exon'][e4][3]
					print('located at ' + str(dis2) + ' bp of closest flanking exon on the left side (' + gff['exon'][e4][8] + ') and ' + str(dis1) + 'bp of the closest flanking exon on the right side (' + gff['exon'][b4][8] +')' )

	else:
		print(pos, 'is intergenic')
		#print('is intergenic')
		#print(gff['gene'][b])
		#print(gff['gene'][e])
		#dis1 = gff['gene'][b][2] - pos
		#dis2 = pos - gff['gene'][e][3]
		#print('located at ' + str(dis2) + ' bp of closest flanking gene on the left side (' + gff['gene'][e][6] + ') and ' + str(dis1) + 'bp of the closest flanking exon on the right side (' + gff['gene'][b][6] +')' )

def check_mutation3(row, arg):
	pos = int(row['POS'])
	print(pos)
	gff = arg['gff']
	b,e = find_row(gff['gene'], pos)
	b15,e15 = find_row(gff['ncRNA_gene'], pos)
	b6,e6 = find_row(gff['ncRNA'], pos)
	b7,e7 = find_row(gff['lnc_RNA'], pos)
	b8,e8 = find_row(gff['snoRNA'], pos)
	b9,e9 = find_row(gff['snRNA'], pos)
	b10,e10 = find_row(gff['rRNA'], pos)
	b11,e11 = find_row(gff['pseudogene'], pos)
	b12,e12 = find_row(gff['pseudogenic_transcript'], pos)
	b13,e13 = find_row(gff['pre_miRNA'], pos)
	b14,e14 = find_row(gff['SRP_RNA'], pos)
	if gff['gene'][b:e+1]:
		print('is genic')
		for gene in gff['gene'][b:e+1]:
			b1,e1 = find_row(gff['mRNA'], pos)
			if gff['mRNA'][b1:e1+1]:
				print('is mRNA')
				print(gff['mRNA'][b1:e1+1])
				for mRNA in gff['mRNA'][b1:e1+1]:
					
					print(mRNA)
					b2,e2 = find_row(gff['exon'], pos)
					if gff['exon'][b2:e2+1]:
						print('is exon')
						b3,e3 = find_row(gff['CDS'], pos)
						print(gff['CDS'][b3:e3+1])
						if gff['CDS'][b3:e3+1]:
							print('is CDS')

							for cds in gff['CDS'][b3:e3+1]:
								
								beg, end = codon_coords(pos,cds[2],cds[3],cds[4],cds[5]) #target pos, cds start, cds end, cds strand, cds phase
								if beg < cds[2]:
									print("  Codon is truncated at exon/intron boundary at: ",cds[2])
									print("  Possible splice site mutation. Check adjacent exon for a possible aa substitution")
								elif end > cds[3]:
									print("  Codon is truncated at exon/intron boundary at: ",cds[3])
									print("  Possible splice site mutation. Check adjacent exon for a possible aa substitution")
								else:
									print(' Codon fully contained exon')
									codon = arg['ref'].seq[beg-1:end]
									before = codon[:pos-beg]+row['REF']+codon[pos-beg+1:]
									after = codon[:pos-beg]+row['ALT']+codon[pos-beg+1:]
									print(before, after)
									find_effect(before,after,cds[4])
						else:
							b4,e4 = find_row(gff['five_prime_UTR'], pos)
							b5,e5 = find_row(gff['three_prime_UTR'], pos) 
							if gff['five_prime_UTR'][b4:e4+1]:
								print('is 5 UTR')
								for five in gff['five_prime_UTR'][b4:e4+1]:

									if five[4] == '+':
										distfive1 = five[3] - pos
										print(' Dist:', distfive1)
									elif five[4] == '-':
										distfive1 = pos - five[2]
										print(' Dist:', distfive1)
							if gff['three_prime_UTR'][b5:e5+1]:
								print(' is 3 UTR')
					else:
						print('is intron')
						b3,e3 = find_row(gff['CDS'], pos)
						dis1 = gff['exon'][b2][2] - pos
						dis2 = pos - gff['exon'][e2][3]
						print('located at ' + str(dis2) + ' bp of closest flanking exon on the left side (' + gff['exon'][e2][6] + ') and ' + str(dis1) + 'bp of the closest flanking exon on the right side (' + gff['exon'][b2][6] +')' )
	elif gff['ncRNA_gene'][b15:e15+1]:
		print('is non coding gene')
	elif gff['ncRNA'][b6:e6+1]:
		print('is non coding RNA')
	elif gff['lnc_RNA'][b7:e7+1]:
		print('is a long non coding RNA')
	elif gff['snoRNA'][b8:e8+1]:
		print('is sno RNA')
	elif gff['snRNA'][b9:e9+1]:
		print('is sn RNA')
	elif gff['rRNA'][b10:e10+1]:
		print('is ribosomic RNA')
	elif gff['pseudogene'][b11:e11+1]:
		print('is a pseudogene')
	elif gff['pseudogenic_transcript'][b12:e12+1]:
		print('is pseudogenic transcript')
	elif gff['pre_miRNA'][b13:e13+1]:
		print('is a precursor of miRNA')
	elif gff['SRP_RNA'][b14:e14+1]:
		print('is a SRP_RNA')

	else:
		print('is intergenic')
		print(gff['gene'][b])
		print(gff['gene'][e])
		dis1 = gff['gene'][b][2] - pos
		dis2 = pos - gff['gene'][e][3]
		print('located at ' + str(dis2) + ' bp of closest flanking gene on the left side (' + gff['gene'][e][6] + ') and ' + str(dis1) + 'bp of the closest flanking exon on the right side (' + gff['gene'][b][6] +')' )


def codon_coords(pos, start, end, strand, phase):
   if strand == '+':
       rest = (pos - start - phase) % 3
       #print("Hebra: + Posicion: "+str(pos))
       #print("Start: "+str(pos-resto), "End: "+str(pos-resto+2))
       return (pos-rest, pos-rest+2)
   elif strand == '-' :
       rest = (end - pos - phase) % 3
       #print("Hebra: - Posicion: "+str(pos))
       #print("Start: "+str(pos+resto-2), "End: "+str(pos+resto))
       return (pos+rest-2, pos+rest)

def get_nearest(arg, pos, type_):
	d = arg['gff'].df
	df_ch = d[d['seq_id'] == arg['gff_chrom']]
	df_ty = df_ch[df_ch['type'] == type_]

	#Nearest object at left of the mutation
	ty_ends = df_ty[(df_ty['end'] <= pos)]
	left_ty = ty_ends.iloc[(ty_ends['end']-pos).abs().argsort().iloc[:1]]
	left = left_ty['end'].loc[left_ty.index[0]]
	lab_left = left_ty['attributes'].loc[left_ty.index[0]].split(';')[0]

	#Nearest object at right to the mutation
	ty_starts = df_ty[(df_ty['start'] >= pos)]# & (exons['strand'] == '+')]
	right_ty = ty_starts.iloc[(ty_starts['start']-pos).abs().argsort().iloc[:1]]
	right = right_ty['start'].loc[right_ty.index[0]]
	lab_right = right_ty['attributes'].loc[right_ty.index[0]].split(';')[0]

	return [pos-left,lab_left,right-pos,lab_right]

def find_effect(before,after,strand):
	before = Seq(before)
	after = Seq(after)
	if strand == '-':
		before = before.reverse_complement()
		after = after.reverse_complement()
	before = before.translate(table=1)
	after = after.translate(table=1)
	
	if (before == after):
		print("Synonymous",before,after)
	elif (before != after):
		if (after == "*"):
			print("Nonsynonymous - nonsense mutation",before,after)
		elif (before == "*"):
			print("Nonsynonymous - nonstop mutation",before,after)
		else:
			print("Nonsynonymous - missense mutation",before,after)

#first = False
#n_head = spacer.join(header)+'\n'
#write_line(n_head, fsal)

#CHROM	POS		REF		ALT		DPr		DPalt		cRef	cAlt	type	      strand		ID		LEFT	RIGHT	Ldist   Rdist 
# 1    1456654   C		T		1         18	ATG(Met)   CTG(Ser) Nonsynonymous   +          ATG10809 ATG..    ATG..  1500    156