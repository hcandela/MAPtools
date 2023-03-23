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
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime,date
import gzip
#from docopt import docopt

def cmHH(command_line):
	os.system(command_line)

def check_args(__doc__,arg:dict):
	'''Check mbs, qtl and annotate arguments.
	   		·Read input
			·Split data
			·Check output
	'''
	arg['--data'] = arg['--data'].split(',')
	if arg['--input'] == None and arg['pipe'] == True:
		print(__doc__, end='\n', file=sys.stderr)
		sys.exit()
	if arg['--input'] != None:
		inp_f = arg['--input']
		if arg['--input'].split('.')[-1] == 'gz':
			try:
				f = gzip.open(inp_f, 'rt')
			except FileNotFoundError:
				print('Error: The input file {} does not exist'.format(inp_f), file=sys.stderr)
				sys.exit()
		else:
			try:
				f = open(inp_f, 'r')
			except FileNotFoundError:
				print('Error: The input file {} does not exist'.format(inp_f), file=sys.stderr)
				sys.exit()
		arg['inp'] = f
	if arg['--output-type'] in {'csv','txt'}:
		if arg['--output-type'] == 'csv':
			arg['spacer'] = ','
		elif arg['--output-type'] == 'txt':
			arg['spacer'] = '\t'
		arg['--output-type'] = '.'+arg['--output-type']
	else:
		print('Error: select a valid format.', file=sys.stderr)
		sys.exit()
	if arg['--output'] != None:
		wd = os.getcwd()
		outdir_list = arg['--output'].split('/')[:-1]
		outdir = '/'.join(outdir_list) +'/'
		try:
			os.makedirs(wd+'/'+outdir)
		except FileExistsError:
			print('Warning: the output directory already exists', file=sys.stderr)
			pass
		arg['filename'] = arg['--output'].split('/')[-1]
		arg['outdir']=wd+'/'+outdir
	
		arg['--output'] = arg['outdir'] + arg['filename']
		arg['--output'] = check_save_an(arg, arg['filename'])
	#Data split, when Wr or Wd is present
	data = arg['--data']
	data_w = data.copy()
	wt_idx = -1
	wt = False
	if 'Wr' in data_w or 'Wd' in data_w:
		if 'Wr' in data_w:
			wt_idx = data_w.index('Wr')
			wt = 'Wr'
		elif 'Wd' in data_w:
			wt_idx = data_w.index('Wd')
			wt = 'Wd'
		data_w.pop(wt_idx)
	arg['data_w'] = data_w
	arg['wt'] = wt
	
	arg['lim'] = 1e-90
	#Call to specific checking function for each command
	if 'mbs' in arg.keys():
		arg = check_mbs_args(arg)
	if 'qtl' in arg.keys():
		arg = check_qtl_args(arg)
	if 'annotate' in arg.keys():
		arg = check_annotate_args(arg)
	return arg

		
def check_save_an(arg:dict, file_name:str):
	"""Check if file exists and make a copy with different name
	"""
	typ = arg['--output-type']
	if os.path.isfile(arg['--output']):
		expand = 1
		while True:
			expand += 1
			nw_file_name = file_name.split(typ)[0] + '(' + str(expand) + ')' + typ
			if os.path.isfile(arg['outdir']+nw_file_name):
				continue
			else:
				file_name = arg['outdir']+ nw_file_name
				return file_name
	else:
		return arg['outdir']+file_name

def check_mbs_args(arg:dict):
	"""Check mbs options
	   		·Also it is used in annotate command
	"""
	data = arg['--data']
	if ('Pr' in data or 'Pd' in data) and arg['wt'] == True:
		#print(wt)
		print('Error: Do not include wild-type(\"Wr\"|\"Wd\") and parental re-sequencing (\"Pd\"|\"Pr\") in the same analysis.', file=sys.stderr)
		sys.exit()
	if arg['--mutant-pool'] not in {'R','D'}:
		print('Error: Select valid mutant pool -m (\"R\"|\"D\")', file=sys.stderr)
		sys.exit()

	if not 'R' in data:
		print('Error: You should include the recessive pool (--data R,X,X)', file=sys.stderr)
		sys.exit()
	if 'Pr' in data and 'Pd' in data:
		print('Error: You should include only one re-sequenced parental in data (--data D,R,Px or --data R,Px)', file=sys.stderr)
		sys.exit()
	if 'Wr' in data and 'Wd' in data:
		print('Error: You should include only one re-sequenced wilt-type in data (--data D,R,Wx or --data R,Wx)', file=sys.stderr)
		sys.exit()
	if arg['--max-depth'] == 'inf':
		arg['--max-depth'] = np.inf
	else:
		arg['--max-depth'] = int(arg['--max-depth'])
	arg['--min-depth'] = int(arg['--min-depth'])
	arg['--min-ratio'] = float(arg['--min-ratio'])/100
	arg['--max-ratio'] = float(arg['--max-ratio'])/100
	arg['--min-error'] = int(arg['--min-error'])
	if arg['--min-depth'] <= 0:
		arg['--min-depth'] = 1
	if arg['--max-depth'] <= arg['--min-depth']:
		print('Error: You should choose a correct interval of depths(--min_depth 20 --max-depth 120)', file=sys.stderr)
		sys.exit()
	if arg['--min-ratio'] < 0:
		arg['--min-ratio'] = 0.0
	if arg['--max-ratio'] <= arg['--min-ratio']:
		print('Error: You should choose a correct interval of frequencies(--min_ratio 15 --max-ratio 85)', file=sys.stderr)
		sys.exit()
	return arg

def check_qtl_args(arg:dict):
	"""Check qtl options"""
	data = arg['--data']
	if len(data) < 2:
		print('Error: You should include a minimum of 2 pools in QTL-seq experiment (--data D,R,Px or --data D,R)', file=sys.stderr)
		sys.exit()
	if 'Pr' in data and 'Pd' in data:
		print('Error: You should include only one parental pool in data (--data D,R,Px)', file=sys.stderr)
		sys.exit()

	arg['--min-depth'] = int(arg['--min-depth'])
	arg['--max-depth'] = int(arg['--max-depth'])
	if arg['--min-depth'] <= 0:
		arg['--min-depth'] = 1
	if arg['--max-depth'] <= arg['--min-depth']:
		print('Error: You should choose a correct interval of depths(--min_depth 20 --max-depth 120)', file=sys.stderr)
		sys.exit()
	return arg

def LF(inp:int):
	"""Calculates factorial as given number, as the summatory of the logarithms (n,1)
	returns a float"""
	result = 0
	if inp == 0:
		return 0
	else:
		for i in range(inp, 1, -1):
			if i >= 1:
				result += log(i)
			else:
				break
	return result

def LogFisher(a:int,b:int,c:int,d:int):
	'''
	Calculates fisher test for p-value calculation.
		·Input: allele depths, int
		·Output: log(fisher_test), float
	'''
	logfisher = (LF(a+b) + LF(c+d) + LF(a+c) + LF(b+d) - LF(a) - LF(b) - LF(c) - LF(d) - LF(a+b+c+d)) / log(10)
	return logfisher

def Gstatic(a:int,b:int,c:int,d:int):
	"""
	Calculates the G statistic for qtl analysis.
		·Input: allele depths, int
		·Output: g, float
	"""
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
	"""Calculates p-value 
	"""
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
	data = arg['--data']
	inf_s = set(data)
	if len(inf_s) == 1 and arg['--ref-genotype'] == 'miss':
		a, b = inp[0], inp[1]
		if a != 0 or b != 0:
			ratio3 = max(a,b)/(a+b)
			boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
			return [ratio3, boost]
	elif (len(inf_s) == 1 and arg['--ref-genotype'] != 'miss') or (len(inf_s) == 2 and 'D' not in inf_s):
		a, b = inp[0], inp[1]
		if a != 0 or b != 0:
			ratio1 = b/(a+b)
			boost = 1/(arg['lim'] + abs(1 - 1/max(ratio1, 1-ratio1)))
			return [ratio1, boost]
	elif 'D' in inf_s and 'R' in inf_s:
		a, b, c, d = inp[0], inp[1], inp[2], inp[3]
		ratio3 = max(c,d)/(c+d)
		resultado = LogFisher(a,b,c,d)
		pva = pvalor(a,b,c,d)
		pva10 = (log(pva)/log(10))
		if arg['--ref-genotype'] == 'miss' and ('Pr' not in inf_s and 'Pd' not in inf_s):
			boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
			return [ratio3,resultado,boost,pva,pva10]
		elif arg['--ref-genotype'] != 'miss' or 'Pr' in inf_s or 'Pd' in inf_s:
			ratio1 = b/(a+b)
			ratio2 = d/(c+d)
			boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
			return [ratio1,ratio2,ratio3,resultado,boost,pva,pva10]

def qtl_calc(inp, arg):
	data = arg['--data']
	inf_s = set(data)
	no_ref = arg['--ref-genotype']
	min_dp = arg['--min-depth']
	max_dp = arg['--max-depth']
	a,b,c,d = inp[0], inp[1], inp[2], inp[3]
	dom = a + b
	rec = c + d
	dp_total = dom + rec
	if dp_total <= max_dp and dp_total >= min_dp and rec > 0 and dom > 0:
		if b/(a+b) == d/(c+d):
			return None
		v1 = (c/rec, d/rec)
		v2 = (a/dom, b/dom)
		ed = distance.euclidean(v1,v2)
		g = Gstatic(a,b,c,d)
		pva = pvalor(a,b,c,d)
		pva10 = (log(pva)/log(10))
		if no_ref == 'miss':
			return [ed,g,pva,pva10]
		else:
			ratio1 = b/(a+b)
			ratio2 = d/(c+d)
			delta = ratio1 - ratio2
			return [ratio1,ratio2,delta,ed,g,pva,pva10]

def choose_header(arg):
	data = arg['--data']
	inf_s = set(data)
	if 'mbs' in arg.keys():
		if len(inf_s) == 1 or 'D' not in inf_s:
			if arg['--ref-genotype'] == 'miss' and 'Pr' not in inf_s and 'Pd' not in inf_s:
				header = ['#CHROM','POS','DOM','REC','DPdom_1','DPrec_1','MAX_SNPidx2','BOOST']
			else:
				header = ['#CHROM','POS','DOM','REC','DPdom_1','DPrec_1','SNPidx1','BOOST']
		elif ('D' in inf_s and 'R' in inf_s) and arg['--ref-genotype'] != 'miss' or 'Pr' in inf_s or 'Pd' in inf_s:
			header = ['#CHROM','POS','DOM','REC','DPdom_1','DPrec_1','DPdom_2','DPrec_2','SNPidx1','SNPidx2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE']
		elif 'D' in inf_s and 'R' in inf_s and arg['--ref-genotype'] == 'miss' and 'Pd' not in inf_s and 'Pr' not in inf_s:
			header = ['#CHROM','POS','DOM','REC','DPdom_1','DPrec_1','DPdom_2','DPrec_2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE']
	elif 'qtl' in arg.keys():
		if arg['--ref-genotype'] == 'miss':
			header = ['#CHROM','POS','DOM','REC','DPdom_1','DPrec_1','DPdom_2','DPrec_2','ED','G','PVALUE','log10PVALUE']
		if arg['--ref-genotype'] != 'miss':
			header = ['#CHROM','POS','DOM','REC','DPdom_1','DPrec_1','DPdom_2','DPrec_2','SNPidx1','SNPidx2','DELTA','ED','G','PVALUE','log10PVALUE']
	arg['header'] = header


def new_line(fsal, arg, first, fields, al_count, calcs):
	res = fields + al_count + calcs
	spacer = arg['spacer']
	header = arg['header']
	if first:	
		first = False
		n_head = spacer.join(header)+'\n'
		write_line(n_head, fsal)

	n_line = spacer.join(str(field) for field in res) + '\n'
	write_line(n_line, fsal)
	return first


def write_line(n_line, fsal):
	if fsal == False:
		print(n_line, end='', file=sys.stdout)
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
	genotype = dict()
	if alt == '.':
		#if the alternative allele is '.' indicates that there is no variation in the aligned sequences, so that position is discarded
		return 0,0,0
	else:
		form = z[8].split(':')
		GT_index = form.index('GT')
		AD_index = form.index('AD')
		data = arg['--data']
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
			GTrec = data1[GT_index].split('/')
			genotype['R'] = GTrec
			if '.' in set(GTrec):
				return 0,0,0
			ref = z[3]
			alt = z[4].split(',')			#alt could be more than one nt
			pool[ref] = int(AD[0])			#generete a dict for the pool {'nt':AD1, 'nt':AD2 }
			for i in range(len(alt)):
				pool[alt[i]] = int(AD[i+1])
			res[data[0]] = pool				#Save the pool dict naming the pool {'R':{'nt1':AD1, 'nt2':AD2 }}
		elif len(inf_s) >= 2:
			if 'D' in inf_s and 'R' in inf_s:
				data1 = z[inf['D']].split(':')	#save info from each pool
				data2 = z[inf['R']].split(':')
				GT1 = data1[GT_index].split('/')
				GT2 = data2[GT_index].split('/')
				gt = set(GT1+GT2)				#compare the genotypes
				if '.' in gt:#if len(gt) == 1 or '.' in gt: #if all the genotypes are identical the line is discarded
					return 0,0,0
				for p,c in inf.items():		#for each pool and column
					ref = z[3]
					alt = z[4].split(',')
					AD = z[c].split(':')[AD_index].split(',')
					GT = z[c].split(':')[GT_index].split('/')
					pool = {ref:int(AD[0])}
					for i in range(len(alt)):
						pool[alt[i]] = int(AD[i+1])
					res[p] = pool	#{'R':{'nt1':AD1, 'nt2':AD2 }, 'D':{'nt1':AD1, 'nt2':AD2}...}
					genotype[p] = GT
			elif 'D' not in inf_s:
				ref = z[3]
				alt = z[4].split(',')
				data1 = z[inf['R']].split(':')
				GT1 = data1[GT_index].split('/')
				gt = set(GT1)
				if '.' in gt:
					return 0,0,0
				for p,c in inf.items():
					AD = z[c].split(':')[AD_index].split(',')
					GT = z[c].split(':')[GT_index].split('/')
					pool = {ref:int(AD[0])}
					for i in range(len(alt)):
						pool[alt[i]] = int(AD[i+1])
					res[p] = pool
					genotype[p] = GT
					
		return [z[0], z[1], z[3]],res,genotype #returns chromosome, position, reference allele, and the data for each bam)

def filter_mbs(arg, al_count, p_al_count, genotype):
	flag = [True]
	inf_s = set(arg['--data'])
	min_dp = arg['--min-depth']
	max_dp = arg['--max-depth']
	if 'R' in inf_s and 'D' in inf_s:
		REF,ALT,a,b,c,d = al_count[0], al_count[1], al_count[2], al_count[3], al_count[4], al_count[5]
		dom = a + b
		rec = c + d
		if dom <= min_dp or dom >= max_dp or rec <= min_dp or rec >= max_dp: #TODO
			flag.append(False)
		flag.append(genotype_filter(arg, genotype))
		if arg['--het-filter']:
			#flag.append(het_filter(arg,a,b))
			flag.append(het_filter2(arg,genotype))
		if arg['--EMS']:
			flag.append(filter_EMS(arg, REF, ALT))
		if arg['--parental-filter'] and ('Pr' in inf_s or 'Pd' in inf_s or 'Wr' in inf_s or 'Wd' in inf_s):
			#e,f = p_al_count[0], p_al_count[1]
			#flag.append(isogenic_filter(arg,c,d,e,f))
			flag.append(parental_filter(arg,genotype))

	elif 'D' not in inf_s:
		REF,ALT,c,d = al_count[0], al_count[1], al_count[2], al_count[3]
		if c+d <= min_dp or c+d >= max_dp:
			flag.append(False)
		flag.append(genotype_filter(arg, genotype))
		if arg['--EMS'] and (arg['--ref-genotype'] != 'miss' or 'Pr' in inf_s or 'Pd' in inf_s):
			flag.append(filter_EMS(arg, REF, ALT))
		if arg['--parental-filter'] and ('Pr' in inf_s or 'Pd' in inf_s or 'Wr' in inf_s or 'Wd' in inf_s):
			#e,f = p_al_count[0], p_al_count[1]
			#flag.append(isogenic_filter(arg,c,d,e,f))
			flag.append(parental_filter(arg,genotype))
	if False in flag:
		return False
	else:
		return True

def genotype_filter(arg, genotype):
	inf_s = set(arg['--data'])
	GT_rec = set(genotype['R'])
	flag = list()
	if 'D' in inf_s and 'R' in inf_s:
		GT_dom = set(genotype['D'])
		rest = inf_s - {'D','R'}
		if len(rest) > 0:
			flag.append(False if GT_dom == GT_rec == {'0'} else True)
			flag.append(False if GT_dom == GT_rec == {'1'} else True)
		elif len(rest) == 0:
			flag.append(False if GT_dom == GT_rec == {'1'} else True)
	if False in flag:
		return False
	else:
		return True


def het_filter(arg,a,b):
	min_dp = arg['--min-depth']
	max_dp = arg['--max-depth']
	min_ratio = arg['--min-ratio']
	max_ratio = arg['--max-ratio']
	dom = a + b
	if dom <= max_dp and dom >= min_dp and (a/dom) <= max_ratio and (a/dom) >= min_ratio:
		return True
	else:
		return False

def het_filter2(arg,genotype):
	GT_dom = set(genotype['D'])
	if GT_dom != {'0','1'}:
		return False
	else:
		return True


def isogenic_filter(arg,c,d,e,f):
	min_dp = arg['--min-depth']
	min_ratio = arg['--min-ratio']
	max_ratio = arg['--max-ratio']
	min_error = arg['--min-error']
	if c + d >= min_dp and e + f >= min_dp and  min_ratio*c < min_dp < max_ratio*c and min_ratio*e < min_dp < max_ratio*e:
		return False
	elif c == c+d and e >= min_error:
		return False
	else: 
		return True
	
def parental_filter(arg,genotype):
	inf_s = set(arg['--data'])
	GT_rec = set(genotype['R'])
	flag = list()
	if 'D' in inf_s and 'R' in inf_s: #Si tenemos 3 pools
		GT_dom = set(genotype['D'])
		rest = inf_s - {'D','R'}
		if len(rest) > 0:
			p = rest.pop()
			GT_par = set(genotype[p])
			flag.append(False if GT_par == GT_dom == GT_rec else True)
	elif 'D' not in inf_s: #Si tenemos 2 pools
		rest = inf_s - {'R'}
		p = rest.pop()
		GT_par = set(genotype[p])
	#Si tenemos 2 o 3 pools y uno es parental o wild-type
	if 'Pr' in inf_s or 'Wr' in inf_s:
		flag.append(False if GT_par == {'0','1'} else True)
		if arg['--mutant-pool'] == 'D':
			flag.append(False if '1' in GT_par else True)
	#elif 'Wr' in inf_s or 'Wd' in inf_s:
	#	flag.append(False if GT_par == {'0','1'} or GT_par == {'1'} else True)
	elif ('Pd' in inf_s or 'Wd' in inf_s) and arg['--mutant-pool'] == 'R':
		flag.append(False if '1' in GT_par else True)
		#TODO-flag.append(False if GT_par == GT_rec else True)
	elif ('Pd' in inf_s or 'Wd' in inf_s) and arg['--mutant-pool'] == 'D':
		flag.append(False if GT_par == {'0'} else True)
	if False in flag:
		return False
	else:
		return True

def outcross_filter(arg,c,d,e,f):
	#if ('Wr' in arg['--data'] or 'Pr' in arg['--data']) and arg['--mutant-pool'] == 'R': #TODO
	if f >= arg['--min-error']:
		return False
	else:
		return True

def normalize(pools, REF, arg, r_min=0.03):
	data = arg['data_w']
	wt = arg['wt']
	inf_s = set(data)
	ref = arg['--ref-genotype']
	rec = pools['R']
	alle = [a for a,v in rec.items()]
	#REF_idx = alle.index(REF)
	arg['reorder'] = 0
	al_count = 0
	wt_l = 0
	if len(alle) > 2:
		return 0,0
	if len(inf_s) == 1:
		if ref == 'D' or ref == 'miss':
			a = rec[alle[0]]
			b = rec[alle[1]]
			al_count = [REF, alle[1], a, b]
		elif ref == 'R':
			b = rec[alle[0]]
			a = rec[alle[1]]
			al_count = [alle[1],REF,a, b]
			arg['reorder'] = 1
	elif len(inf_s) == 2:
		if 'R' in inf_s and 'D' in inf_s and ref == 'miss':
			dom = pools['D']
			a = dom[alle[0]]
			c = rec[alle[0]]
			b = dom[alle[1]]
			d =	rec[alle[1]]
			al_count = [alle[0], alle[1], a, b, c, d]
		elif inf_s == {'R','D'} and ref == 'D':
			dom = pools['D']
			a = dom[alle[0]]
			c = rec[alle[0]]
			b = dom[alle[1]]
			d = rec[alle[1]]
			al_count = [REF, alle[1], a, b, c, d]
		elif inf_s == {'R','D'} and ref == 'R':
			dom = pools['D']
			b = dom[alle[0]]
			d = rec[alle[0]]
			a = dom[alle[1]]
			c = rec[alle[1]]
			al_count = [alle[1],REF, a, b, c, d]
			arg['reorder'] = 1
		elif inf_s == {'R', 'Pd'}:
			p_dom = pools['Pd']
			p_al = sorted(p_dom, key=lambda key: p_dom[key], reverse=True) #ordered from higher to lower
			p_al2 = [key for key,arg in p_dom.items()]
			arg['reorder'] = 1 if p_al != p_al2 else 0
			if p_dom[p_al[0]] > 0 and p_dom[p_al[1]]/(p_dom[p_al[0]] + p_dom[p_al[1]]) < r_min:# or len(p_al) > 2:
				a = rec[p_al[0]]
				b = rec[p_al[1]]
				e = p_dom[p_al[0]]
				f = p_dom[p_al[1]]
				al_count = [p_al[0], p_al[1], a, b]
				wt_l = [e,f]
		elif inf_s == {'R', 'Pr'}:	#TODO Possibly delete this case
			p_rec = pools['Pr']
			p_al = sorted(p_rec, key=lambda key: p_rec[key], reverse=True)
			p_al2 = [key for key,arg in p_rec.items()]
			arg['reorder'] = 1 if p_al != p_al2 else 0
			if p_rec[p_al[0]] > 0 and p_rec[p_al[1]]/(p_rec[p_al[0]] + p_rec[p_al[1]]) < r_min:# or len(p_al) > 2:
				a = rec[p_al[1]]
				b = rec[p_al[0]]
				e = p_rec[p_al[1]]
				f = p_rec[p_al[0]]
				al_count = [p_al[1], p_al[0], a, b]
				wt_l = [e, f]
			else:
				al_count, wt_l = 0,0
	elif len(inf_s) == 3:
		if 'Pr' in inf_s:
			parental = pools['Pr']
		elif 'Pd' in inf_s:
			parental = pools['Pd']
		p_al = sorted(parental, key=lambda key: parental[key], reverse=True)
		dom = pools['D']
		r_min = 0.03
		if parental[p_al[0]] > 0 and parental[p_al[1]]/(parental[p_al[0]] + parental[p_al[1]]) < r_min:
			p_al2 = [key for key,arg in parental.items()]
			arg['reorder'] = 1 if p_al != p_al2 else 0
			if 'Pd' in inf_s:
				a = dom[p_al[0]]
				b = dom[p_al[1]]
				c = rec[p_al[0]]
				d = rec[p_al[1]]
				e = parental[p_al[0]]
				f = parental[p_al[1]]
				al_count = [p_al[0], p_al[1], a, b, c, d]
				wt_l = [e, f] 
			elif 'Pr' in inf_s:
				a = dom[p_al[1]]
				b = dom[p_al[0]]
				c = rec[p_al[1]]
				d = rec[p_al[0]]
				e = parental[p_al[1]]
				f = parental[p_al[0]]
				al_count = [p_al[1], p_al[0], a, b, c, d]
				wt_l = [e, f]
		else:
			al_count, wt_l = 0, 0
	if wt:
		if arg['reorder'] == 1:
			wt_l = [pools[wt][alle[1]], pools[wt][alle[0]]]
		elif arg['reorder'] == 0:
			wt_l = [pools[wt][alle[0]], pools[wt][alle[1]]]
	#elif not wt and ('Pr' not in inf_s or 'Pd' not in inf_s):
	#	wt_l = 0
	return al_count, wt_l
	

def check_annotate_args(arg):
	arg['--version'] = True
	wd = os.getcwd()
	if arg['--gff'] == None:
		print('Error: Select a gff file (--gff FILE)', file=sys.stderr)
		sys.exit()
	try:
		if arg['--gff'].split('.')[-1] == 'gz':
			gff = gzip.open(wd+'/'+arg['--gff'],'rt')
		else:
			gff = open(wd+'/'+arg['--gff'],'r')
		arg['--gff'] = wd+'/'+arg['--gff']
		arg['--gff'] = gff
	except FileNotFoundError:
		print('The gff file does not exist.', file=sys.stderr)
		sys.exit()

	if arg['--fasta-reference'] == None:
		print('Error: Select a fasta reference file (--fasta-reference FILE)', file=sys.stderr)
		sys.exit()
	try:
		if arg['--fasta-reference'].split('.')[-1] == 'gz':
			ref = gzip.open(wd+'/'+arg['--fasta-reference'],'rt')
		else:
			ref = open(wd+'/'+arg['--fasta-reference'],'r')
		arg['--fasta-reference'] = wd+'/'+arg['--fasta-reference']
		ref.close()
	except FileNotFoundError:
		print('The reference file does not exist.', file=sys.stderr)
		sys.exit()
	#arg['chromosomes'] = list()
	if arg['--region'] == None:
		print('Error: Select a region (-R chromName:startPos-endPos)', file=sys.stderr)
		sys.exit()
	else:
		if re.match('\S+:\d+-\d+', arg['--region']):
			Tchrom,Treg = arg['--region'].split(':')
			Treg = Treg.split('-')
			arg['--region'] = [Tchrom,int(Treg[0]),int(Treg[1])]
		else:
			print('Error: Enter region correctly (-R chromName:startPos-endPos)', file=sys.stderr)
			sys.exit()
	if 'R' in arg['--data'] and 'D' in arg['--data']:
		header=['#CHROM','POS','DOM','REC', 'DPdom_1','DPrec_1','DPdom_2','DPrec_2','TYPE','ID','PARENT','STRAND',\
			'PHASE','CODON_change','AA_change','INFO']
	else:
		header=['#CHROM','POS','DOM','REC', 'DPdom_1','DPrec_1','TYPE','ID','PARENT','STRAND',\
			'PHASE','CODON_change','AA_change','INFO']
	arg['header'] = header
	if not 'R' in arg['--data']:
		print('Error: You should include the recessive pool (--data R,X)', file=sys.stderr)
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
	#with open(gff_path, 'r') as handle:
	for line in gff_path:
		if line.startswith('#'):
			continue
		else:
			line = line.rstrip()
			line = line.rstrip(';')
			data = line.split('\t')
			seqid = data[0]
			#TODO Complete condition to cath only structures in our region
			if seqid == arg['--region'][0]: #Load only the chromosome selected
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
					if type_ == 'CDS' and phase == '.':
						print('##Warning: Phase has to be an integer in CDS {}:{}-{} and was discarded.\n'.format(seqid,start,end), file=sys.stderr)
					else:			
						gff[type_].append([seqid, type_, int(start), int(end), strand, phase, dict_att['ID'], dict_att['Parent'], dict_att['Name']])

	gff['gene'] = gff['gene'] + gff['ncRNA_gene'] + gff['pseudogene']	#all of these types are considered as genes
	gff['gene'].sort(key=lambda row:row[2]) #ordered by pos
	gff['mRNA'].sort(key=lambda row:(row[6], row[2])) #ordered by id and pos
	gff['exon'].sort(key=lambda row:(row[7], row[2])) #ordered by parent and pos
	gff['CDS'].sort(key=lambda row:(row[7], row[2])) #ordered by parent and pos
	gff['five_prime_UTR'].sort(key=lambda row:(row[7], row[2])) #ordered by parent and pos
	gff['three_prime_UTR'].sort(key=lambda row:(row[7], row[2])) #ordered by parent and pos
	ordered = ['gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']
	for type_ in gff: #the rest of types are ordered only by pos
		if type_ not in ordered:
			gff[type_].sort(key=lambda row:row[2])
	arg['gff'] = gff

def codon_coords(pos, start, end, strand, phase):
	if strand == '+':
		rest = (pos - start - phase) % 3
		return (pos-rest, pos-rest+2)
	elif strand == '-' :
		rest = (end - pos - phase) % 3
		return (pos+rest-2, pos+rest)
				
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
	if arg['--fasta-reference'].split('.')[-1] == 'gz':
		with gzip.open(arg['--fasta-reference'], 'rt') as handle:
			for sequ in SeqIO.parse(handle, 'fasta'):
				if sequ.id == chrom:
					arg['ref'] = sequ
					flag = False
					break
			if flag == True:
				print('Error: Chromosome not found in FASTA file.', file=sys.stderr)
				sys.exit()
	else:
		with open(arg['--fasta-reference'], 'r') as handle:
			for sequ in SeqIO.parse(handle, 'fasta'):
				if sequ.id == chrom:
					arg['ref'] = sequ
					flag = False
					break
			if flag == True:
				print('Chromosome not found in FASTA file.', file=sys.stderr)
				sys.exit()


def find_effect(before,after,strand):
	before = Seq(before)
	after = Seq(after)
	if strand == '-':
		before = before.reverse_complement()
		after = after.reverse_complement()
	before = before.translate(table=1)
	after = after.translate(table=1)
	if (before == after):
		change = 'Synonymous:.'
	elif (before != after):
		if (after == '*'):
			change = 'Non_synonymous:nonsense'
		elif (before == '*'):
			change = 'Nonsynonymous:nonstop'
		else:
			change = 'Nonsynonymous:missense'
	return [before, after,change]

def filter_EMS(arg, REF, ALT):
	if arg['--mutant-pool'] == 'R':
		if (REF == 'G' and ALT == 'A') or (REF == 'C' and ALT == 'T'):
			return True
		else:
			return False
	elif arg['--mutant-pool'] == 'D':
		if (REF == 'A' and ALT == 'G') or (REF == 'T' and ALT == 'C'):
			return True
		else:
			return False


def create_df(arg):
	data = arg['--data']
	inf_s = set(data)
	if len(inf_s) == 1 or 'D' not in inf_s:
		if arg['--ref-genotype'] == 'miss':
			header = ['#CHROM','POS','DOM','REC','DPdom_1','DPrec_1','MAX_SNPidx2','BOOST','reorder']
		else:
			header = ['#CHROM','POS','DOM','REC','DPdom_1','DPrec_1','SNPidx1','BOOST','reorder']
	elif ('D' in inf_s and 'R' in inf_s) and arg['--ref-genotype'] != 'miss' or 'Pr' in inf_s or 'Pd' in inf_s:
		header = ['#CHROM','POS','DOM','REC','DPdom_1','DPrec_1','DPdom_2','DPrec_2','SNPidx1','SNPidx2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE','reorder']
	elif 'D' in inf_s and 'R' in inf_s and arg['--ref-genotype'] == 'miss' and 'Pd' not in inf_s and 'Pr' not in inf_s:
		header = ['#CHROM','POS','DOM','REC','DPdom_1','DPrec_1','DPdom_2','DPrec_2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE','reorder']
	arg['header2'] = header
	df = pd.DataFrame(columns=header)
	return df

def new_df_line(df, arg, fields, al_count, calcs):
	data = arg['--data']
	inf_s = set(data)
	reorder = [arg['reorder']]
	res = fields + al_count + calcs + reorder
	n_line = {arg['header2'][i]:res[i] for i in range(len(res))}
	df_nline = pd.DataFrame.from_records([n_line])
	df = pd.concat([df, df_nline])
	return df

def find_new_ATG(arg, pos, five, row):
	if five[4] == '+':
		c1 = arg['ref'].seq[pos-3:pos-1]+row['REC']
		c2 = arg['ref'].seq[pos-2:pos-1]+row['REC']+arg['ref'].seq[pos:pos+1]
		c3 = row['REC']+arg['ref'].seq[pos:pos+2]
		if c1 == 'ATG' or c2 == 'ATG' or c3 == 'ATG':
			if c1 == 'ATG':
				lag = -2
			elif c2 == 'ATG':
				lag = -1
			elif c3 == 'ATG':
				lag = 0
			dist_to_cds = five[3] - (pos + lag) + 1
			if (dist_to_cds)%3 == 0:
				codon = (Seq('ATG') + arg['ref'].seq[five[3]-dist_to_cds+3:five[3]]).translate(table=1)
				if '*' in codon:
					return 'new_ATG:in_frame:truncated_protein'
				else:
					return 'new_ATG:in_frame:elongated_protein'
			else:
				return 'new_ATG:out_of_frame:.'
		else:
			return '.:.:.'
	if five[4] == '-':
		c1 = arg['ref'].seq[pos-3:pos-1]+row['REC']
		c2 = arg['ref'].seq[pos-2]+row['REC']+arg['ref'].seq[pos:pos+1]
		c3 = row['REC']+arg['ref'].seq[pos:pos+2]
		if c1 == 'CAT' or c2 == 'CAT' or c3 == 'CAT':
			if c1 == 'CAT':
				lag = 0
			elif c2 == 'CAT':
				lag = 1
			elif c3 == 'CAT':
				lag = 2
			dist_to_cds = pos + lag - five[2] + 1
			if (dist_to_cds)%3 == 0:
				s_UTR = (arg['ref'][five[2]-1:pos+lag-3]+Seq('CAT')).reverse_complement()
				prot = s_UTR.translate(table=1)
				if '*' in prot:
					return 'new_ATG:in_frame:truncated_protein'
				else:
					return 'new_ATG:in_frame:elongated_protein'
			else:
				return 'new_ATG:out_of_frame:.'
		else:
			return '.:.:.'

def check_nc_gene(arg,gff,type_, b, e, pos,result, row):
	result['PHASE'],result['CODON_ref'], result['CODON_alt'], result['AA_ref'],result['AA_alt']= '.','.','.','.','.'
	if gff[type_][b:e+1]:
		for ncRNA in gff[type_][b:e+1]:
			ncRNA_id = ncRNA[6]
			b1,e1 = find_row_name(gff['exon'], ncRNA_id, 7)
			b2,e2 = find_row(gff['exon'], pos, b1, e1+1)
			if gff['exon'][b2:e2+1]:
				result['ID'],result['PARENT'],result['TYPE'],result['STRAND'] = ncRNA_id,ncRNA[7],type_,ncRNA[4]
				result['INFO']['effect'] = 'non_coding:exonic'
				write_annotate_line(arg, result, row)
				result['INFO'] = dict()
			else:
				result['ID'],result['PARENT'],result['TYPE'],result['STRAND'] = ncRNA_id,ncRNA[7],type_,ncRNA[4]
				result['INFO']['effect'] = 'non_coding:intronic'
				write_annotate_line(arg, result, row)
				result['INFO'] = dict()
			

def check_mutation2(row, arg):
	result = {h:row[h] for h in arg['header2'] if h in arg['header2'] and h in arg['header']}
	result2 = {h:'.' for h in arg['header'] if h not in arg['header2']}
	result2['INFO'] = dict()
	result.update(result2)
	result['CODON_ref'], result['CODON_alt'], result['AA_ref'], result['AA_alt'] = '.','.','.','.'
	gff = arg['gff']
	pos = int(row['POS'])
	b,e = find_row(gff['gene'], pos, 0, len(gff['gene']))
	if gff['gene'][b:e+1]:
		for gene in gff['gene'][b:e+1]:
			gene_id = gene[6]
			b1,e1 = find_row_name(gff['mRNA'], gene_id, 7) #mRNA
			b2,e2 = find_row(gff['mRNA'], pos, b1, e1+1)
			if gff['mRNA'][b2:e2+1]:
				for mRNA in gff['mRNA'][b2:e2+1]:
					mRNA_id = mRNA[6]
					result['STRAND'] = mRNA[4]
					b3,e3 = find_row_name(gff['exon'], mRNA_id, 7)
					b4,e4 = find_row(gff['exon'], pos, b3, e3+1)
					if pos <= mRNA[3] and pos >= mRNA[2]:
						if gff['exon'][b4:e4+1]:
							b5,e5 = find_row_name(gff['CDS'], mRNA_id, 7)
							b6,e6 = find_row(gff['CDS'], pos, b5, e5+1)
							b7,e7 = find_row_name(gff['five_prime_UTR'], mRNA_id, 7)
							b8,e8 = find_row(gff['five_prime_UTR'], pos, b7, e7+1)
							b9,e9 = find_row_name(gff['three_prime_UTR'], mRNA_id, 7)
							b10,e10 = find_row(gff['three_prime_UTR'], pos, b9, e9+1)
							if gff['CDS'][b6:e6+1]:
								for cds in gff['CDS'][b6:e6+1]:
									cds_id = cds[6]
									result['PHASE'],result['TYPE'] = cds[5],'CDS'
									beg, end = codon_coords(pos,cds[2],cds[3],cds[4],int(cds[5])) #target pos, cds start, cds end, cds strand, cds phase
									if gff['CDS'][b6-1][7] == cds[7]: #if prev cds exist
										if pos - cds[2] in {0,1,2}:
											if cds[4] == '+':
												result['INFO']['3_splice_site'] = 'exon-boundary:'+gff['CDS'][e6][6]
											else:
												result['INFO']['5_splice_site'] = 'exon-boundary:'+gff['CDS'][e6][6]

									if gff['CDS'][b6+1][7] == cds[7]: #if next cds exist
										if cds[3] - pos in {0,1,2}:
											if cds[4] == '+':
												result['INFO']['5_splice_site'] = 'exon-boundary:'+gff['CDS'][e6][6]
											else:
												result['INFO']['3_splice_site'] = 'exon-boundary:'+gff['CDS'][e6][6]

									if beg < cds[2]:
										prev_cds = gff['CDS'][e6-1]
										p_beg = prev_cds[2]
										p_end = prev_cds[3]
										seq_cds_p = arg['ref'].seq[p_beg-1:p_end]
										cds = gff['CDS'][e6]
										n_beg = gff['CDS'][e6][2]
										n_end = gff['CDS'][e6][3]
										rest = cds[2] - beg
										codon_ref = seq_cds_p[-rest:] + arg['ref'].seq[n_beg-1: end]
										codon_alt = codon_ref[:pos-beg]+row['REC']+codon_ref[pos-beg+1:]
										aa_before, aa_after,change = find_effect(codon_ref,codon_alt,cds[4])
										result['CODON_ref'],result['CODON_alt'],result['AA_ref'],result['AA_alt'] = codon_ref, codon_alt, aa_before, aa_after
										result['INFO']['effect'] = change
										result['PARENT'] = mRNA_id
										write_annotate_line(arg, result, row)
										result['INFO'] = dict()
										result['PHASE'],result['CODON_ref'], result['CODON_alt'], result['AA_ref'],result['AA_alt']= '.','.','.','.','.'
									elif end > cds[3]:
										next_cds = gff['CDS'][e6+1]
										n_beg = next_cds[2]
										n_end = next_cds[3]
										seq_cds_n = arg['ref'].seq[n_beg-1:n_end]
										p_beg = gff['CDS'][e6][2]
										p_end = gff['CDS'][e6][3]
										rest = end - cds[3]
										codon_ref = arg['ref'].seq[beg-1: p_end] + seq_cds_n[:rest]
										codon_alt = codon_ref[:pos-beg]+row['REC']+codon_ref[pos-beg+1:]
										aa_before, aa_after,change = find_effect(codon_ref,codon_alt,cds[4])
										result['CODON_ref'],result['CODON_alt'],result['AA_ref'],result['AA_alt'] = codon_ref, codon_alt, aa_before, aa_after
										result['PARENT'] = mRNA_id
										result['INFO']['effect'] = change
										write_annotate_line(arg, result, row)
										result['INFO'] = dict()
										result['PHASE'],result['CODON_ref'], result['CODON_alt'], result['AA_ref'],result['AA_alt']= '.','.','.','.','.'
									else:
										codon = arg['ref'].seq[beg-1:end]
										before = codon[:pos-beg]+row['DOM']+codon[pos-beg+1:]
										after = codon[:pos-beg]+row['REC']+codon[pos-beg+1:]
										aa_before, aa_after,change = find_effect(before,after,cds[4])
										result['PARENT'] = mRNA_id
										result['CODON_ref'],result['CODON_alt'],result['AA_ref'],result['AA_alt'] = before, after, aa_before, aa_after
										result['INFO']['effect'] = change
										write_annotate_line(arg, result, row)
										result['INFO'] = dict()
										result['PHASE'],result['CODON_ref'], result['CODON_alt'], result['AA_ref'],result['AA_alt']= '.','.','.','.','.'
							elif gff['five_prime_UTR'][b8:e8+1]:
								#five =  gff['five_prime_UTR'][b8:e8+1][0]
								for five in  gff['five_prime_UTR'][b8:e8+1]:
									result['TYPE'] = 'five_prime_UTR'		
									result['INFO']['effect'] = find_new_ATG(arg, pos, five, row)
									result['PARENT'] = mRNA_id
									write_annotate_line(arg, result, row)
									result['INFO'] = dict()
							elif gff['three_prime_UTR'][b10:e10+1]:
								three =  gff['three_prime_UTR'][b10:e10+1][0]
								result['TYPE'] = 'three_prime_UTR'
								result['PARENT'] = mRNA_id
								write_annotate_line(arg, result, row)
								result['INFO'] = dict()
						else:
							result['TYPE'] = 'intron'
							dis1 = gff['exon'][b4][2] - pos
							dis2 = pos - gff['exon'][e4][3]
							result['PARENT'] = mRNA_id
							result['INFO']['left'] = gff['exon'][e4][8] +':'+ str(dis2)
							result['INFO']['right'] = gff['exon'][b4][8] +':'+ str(dis1)
							if dis2 in {1,2,3}:
								if gff['exon'][e4][4] == '+':
									result['INFO']['5_splice_site'] = 'intron-boundary:'+gff['exon'][e4][8]
								else:
									result['INFO']['3_splice_site'] = 'intron-boundary:'+gff['exon'][e4][8]
							if dis1 in {1,2,3}:
								if gff['exon'][b4][4] == '+':
									result['INFO']['3_splice_site'] = 'intron-boundary:'+gff['exon'][b4][8]
								else:
									result['INFO']['5_splice_site'] = 'intron-boundary:'+gff['exon'][b4][8]
							write_annotate_line(arg, result, row)
							result['INFO'] = dict()
			
			b11,e11 = find_row_name(gff['tRNA'], gene_id, 7) #tRNA
			b12,e12 = find_row(gff['tRNA'], pos, b11, e11+1)
			check_nc_gene(arg,gff,'tRNA', b12, e12, pos,result,row)
			b13,e13 = find_row_name(gff['rRNA'], gene_id, 7) #rRNA
			b14,e14 = find_row(gff['rRNA'], pos, b13, e13+1)
			check_nc_gene(arg,gff,'rRNA', b14, e14, pos,result,row)
			b15,e15 = find_row_name(gff['ncRNA'], gene_id, 7) #ncRNA
			b16,e16 = find_row(gff['ncRNA'], pos, b15, e15+1)
			check_nc_gene(arg,gff,'ncRNA', b16, e16, pos,result,row)
			b17,e17 = find_row_name(gff['lnc_RNA'], gene_id, 7) #lncRNA
			b18,e18 = find_row(gff['lnc_RNA'], pos, b17, e17+1)
			check_nc_gene(arg,gff,'lnc_RNA', b18, e18, pos,result,row)
			b19,e19 = find_row_name(gff['miRNA'], gene_id, 7) #miRNA
			b20,e20 = find_row(gff['miRNA'], pos, b19, e19+1)
			check_nc_gene(arg,gff,'miRNA', b20, e20, pos,result,row)
			b21,e21 = find_row_name(gff['pre_miRNA'], gene_id, 7) #pre_miRNA
			b22,e22 = find_row(gff['pre_miRNA'], pos, b21, e21+1)
			check_nc_gene(arg,gff,'pre_miRNA', b22, e22, pos,result,row)
			b23,e23 = find_row_name(gff['snRNA'], gene_id, 7) #snRNA
			b24,e24 = find_row(gff['snRNA'], pos, b23, e23+1)
			check_nc_gene(arg,gff,'snRNA', b24, e24, pos,result,row)
			b25,e25 = find_row_name(gff['snoRNA'], gene_id, 7) #snoRNA
			b26,e26 = find_row(gff['snoRNA'], pos, b25, e25+1)
			check_nc_gene(arg,gff,'snoRNA', b26, e26, pos,result,row)
			b27,e27 = find_row_name(gff['pseudogenic_transcript'], gene_id, 7) #pseudogenic_transcript
			b28,e28 = find_row(gff['pseudogenic_transcript'], pos, b27, e27+1)
			check_nc_gene(arg,gff,'pseudogenic_transcript', b28, e28, pos,result,row)
			b29,e29 = find_row_name(gff['SRP_RNA'], gene_id, 7) #SRP_RNA
			b30,e30 = find_row(gff['SRP_RNA'], pos, b29, e29+1)
			check_nc_gene(arg,gff,'SRP_RNA', b30, e30, pos,result,row)

	else:
		dis1 = gff['gene'][b][2] - pos
		dis2 = pos - gff['gene'][e][3]
		result['TYPE'] = 'intergenic'
		result['INFO']['left'] = gff['gene'][e][6] +':'+ str(dis2)
		result['INFO']['right'] = gff['gene'][b][6] +':'+ str(dis1)
		write_annotate_line(arg, result, row)
		result['INFO'] = dict()

def write_annotate_header(arg):
	fsal = arg['fsal']
	header = arg['header']
	spacer = arg['spacer']
	n_head = spacer.join(header)+'\n'
	write_line(n_head,fsal)

def write_annotate_line(arg, result, row):
	codonREF,codonALT,AAREF,AAALT =  result['CODON_ref'], result['CODON_alt'], result['AA_ref'], result['AA_alt']
	if codonREF == '.' and codonALT == '.' and AAREF == '.' and AAALT == '.':
		result['CODON_change'] = '.'
		result['AA_change'] = '.'
	else:
		if row['reorder'] == 1:
			result['CODON_change'] = codonALT+' > '+codonREF
			result['AA_change'] = AAALT+' > '+AAREF
		elif row['reorder'] == 0:
			result['CODON_change'] = codonREF+' > '+codonALT
			result['AA_change'] = AAREF +' > '+AAALT

	att_l = [k+'='+v for k,v in result['INFO'].items()]
	att = ';'.join(att_l)
	result['INFO'] = att
	result2 = {field:result[field] for field in arg['header']}
	line = [str(v) for k,v in result2.items()]
	n_line = arg['spacer'].join(line)+'\n'
	write_line(n_line,arg['fsal'])



#CHROM	POS		REF		ALT		DPr		DPalt		cRef	cAlt	type	      strand		ID		LEFT	RIGHT	Ldist   Rdist 
# 1    1456654   C		T		1         18	ATG(Met)   CTG(Ser) Nonsynonymous   +          ATG10809 ATG..    ATG..  1500    156

def read_header(arg:dict,line:str):
	fsal = arg['fsal']
	arg['chromosomes'] = list()
	if line.startswith('##bcftoolsVersion='):
		write_line(line,fsal)
	if line.startswith('##bcftoolsCommand='):
		write_line(line,fsal)
	if line.startswith('##bcftools_callCommand='):
		write_line(line.split(';')[0]+'\n',fsal)
	if line.startswith('##reference='):
		if 'mbs' in arg.keys() or 'qtl' in arg.keys():
			write_line(line,fsal)
	if line.startswith('##contig=<ID='):
		c = line.split('##contig=<ID=')[1].split(',')[0]
		arg['chromosomes'].append(c)
		write_line(line,fsal)
	if line.startswith('#CHROM'):
		bam_list = line.rstrip().split('\t')[9:]
		if len(arg['--data']) > len(bam_list):
			print('Error: The data list does not match with the number of pools used.', file=sys.stderr)
			sys.exit()
		elif len(arg['--data']) <= len(bam_list):
			arg['pools'] = {arg['--data'][i]:bam_list[i] for i in range(len(arg['--data']))}
			pools = '##bamFiles:genotypes='+','.join(':'.join((key,val)) for (key, val) in arg['pools'].items()) + '\n'
			write_line(pools, fsal)

def write_argv(arg:dict,argv:str):
	fsal = arg['fsal']
	version_maptools = '##maptoolsVersion='+v_maptools+'\n'
	version = '##maptools_{}Version='.format(argv[0])+arg['version']+'\n'
	line = '##maptools_{}Command='.format(argv[0])+' '.join(argv) +';'+'Date='+datetime.now().ctime()+'\n'
	write_line(version_maptools,fsal)
	write_line(version,fsal)
	write_line(line, fsal)
	if argv[0] == 'annotate':
		ref_idx = argv.index('--fasta-reference') if '--fasta-reference' in argv else argv.index('-f')
		gff_idx = argv.index('--gff') if '--gff' in argv else argv.index('-g')
		write_line('##reference=file://'+argv[ref_idx+1].split('/')[-1]+'\n',fsal)
		write_line('##gff3=file://'+argv[gff_idx+1].split('/')[-1]+'\n',fsal)
	for field in arg['header']:
		if 'merge' in arg.keys():
			if field == 'DOM' or field == 'REC':
				continue
		if field in variable_descriptions.keys():
			write_line(variable_descriptions[field], fsal)
	if argv[0] == 'merge':
		line ='##maptools_mergeCommandINFO=\"All columns now represent the grouped values for the provided window size.\"\n'
		write_line(line, fsal)

def read_header_merge(arg:dict):
	lines = list()
	with open(arg['--input'], 'r') as handle:
		for line in handle:
			if line.startswith('##maptools_'):
				lines.append(line)
			if line.startswith('##bcftoolsVersion='):
				lines.append(line)
			if line.startswith('##bcftoolsCommand='):
				lines.append(line)
			if line.startswith('##bcftools_callCommand='):
				line = line.split(';')[0]
				lines.append(line)
			if line.startswith('##reference='):
				if 'mbs' in arg.keys() or 'qtl' in arg.keys() or 'merge' in arg.keys():
					lines.append(line)
			if line.startswith('##contig=<ID='):
				c = line.split('##contig=<ID=')[1].split(',')[0]
				#arg['chromosomes'].append(c)
				lines.append(line)
			if line.startswith('##bamFiles'):
				lines.append(line)
			if line.startswith('#CHROM'):
				header = line.split('\t')
				header = [field.rstrip() for field in header]
				arg['header'] = header
				break
	return lines


		