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


def cmHH(command_line):
	os.system(command_line)

def check_args(__doc__,arg:dict):
	'''Check mbs, qtl and annotate arguments.
	   		·Read input
			·Split data
			·Check output
	'''
	arg['--data'] = arg['--data'].split(',')
	if 'mbs' in arg.keys():
		valid_pools = {'D','R','Pr','Pd','Wr','Wd'}
		for d in arg['--data']:
			if d not in valid_pools:
				print(f'Error: \"{d}\" is not a valid sample name.', file=sys.sterr)
				sys.exit()
	if 'qtl' in arg.keys():
		translate = {'H':'D','L':'R','P':'Pd'}
		for d in arg['--data']:
			if d not in translate.keys():
				print(f'Error: \"{d}\" is not a valid sample name.', file=sys.sterr)
				sys.exit()
		arg['--data'] = [translate[i] for i in arg['--data']]
		if arg['--ref-genotype'] not in {'P','miss'}:
			print('Error: Reference genotype\'s tag (-r) should be \'P\' or \'miss\'', file=sys.stderr)
			sys.exit()
		if arg['--ref-genotype'] == 'P':
			arg['--ref-genotype'] = 'D'
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
	arg['--contigs'] = dict()
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
			nw_file_name = file_name.split(typ)[0] + '_' + str(expand) + typ
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
	arg['dp_filter'] = False if arg['--max-depth'] == np.inf and arg['--min-depth'] == 0 else True
	arg['--min-ratio'] = float(arg['--min-ratio'])/100
	arg['--max-ratio'] = float(arg['--max-ratio'])/100
	#arg['--min-error'] = int(arg['--min-error'])
	if arg['--min-depth'] <= 0:
		arg['--min-depth'] = 1
	if arg['--max-depth'] <= arg['--min-depth']:
		print('Error: You should choose a correct interval of depths(--min_depth 20 --max-depth 120)', file=sys.stderr)
		sys.exit()
	if arg['--min-ratio'] < 0:
		arg['--min-ratio'] = 0.0
	if arg['--max-ratio'] <= arg['--min-ratio']:
		print('Error: You should choose a correct interval of frequencies(--min-ratio 15 --max-ratio 85)', file=sys.stderr)
		sys.exit()
	return arg

def check_qtl_args(arg:dict):
	"""Check qtl options"""
	data = arg['--data']
	if 'D' not in data or 'R' not in data:
		print('Error: You should include the 2 pools in QTL-seq experiment (--data H,L)', file=sys.stderr)
		sys.exit()

	if arg['--max-depth'] == 'inf':
		arg['--max-depth'] = np.inf
	else:
		arg['--max-depth'] = int(arg['--max-depth'])
	arg['--min-depth'] = int(arg['--min-depth'])
	arg['dp_filter'] = False if arg['--max-depth'] == np.inf and arg['--min-depth'] == 0 else True
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
		if (a+b) > 0 and (c+d) > 0:
			ratio3 = max(c,d)/(c+d)
			resultado = LogFisher(a,b,c,d)
			pva = pvalor(a,b,c,d)
			pva10 = (log(pva)/log(10))
			if arg['--ref-genotype'] == 'miss' and ('Pr' not in inf_s and 'Pd' not in inf_s):
				boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
				return [ratio3,resultado,boost,pva,pva10]
			elif arg['--ref-genotype'] != 'miss' or 'Pr' in inf_s or 'Pd' in inf_s:
				#if (a+b) > 0 and (c+d) > 0:
				ratio1 = b/(a+b)
				ratio2 = d/(c+d)
				boost = 1/(arg['lim'] + abs(1 - 1/ratio3))
				return [ratio1,ratio2,ratio3,resultado,boost,pva,pva10]

def qtl_calc(inp, arg):
	data = arg['--data']
	inf_s = set(data)
	no_ref = arg['--ref-genotype']
	a,b,c,d = inp[0], inp[1], inp[2], inp[3]
	dom = a + b
	rec = c + d
	if rec > 0 and dom > 0:
		v1 = (c/rec, d/rec)
		v2 = (a/dom, b/dom)
		ed = distance.euclidean(v1,v2)
		g = Gstatic(a,b,c,d)
		pva = pvalor(a,b,c,d)
		pva10 = (log(pva)/log(10))
		ratio1 = b/(a+b)
		ratio2 = d/(c+d)
		delta = ratio1 - ratio2
		if no_ref == 'miss' and 'Pd' not in data:
			return [ratio1, ratio2,abs(delta),ed,g,pva,pva10]
		else:
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
		#if arg['--ref-genotype'] == 'miss' and 'Pr' not in inf_s and 'Pd' not in inf_s:
		#	header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','ED','G','PVALUE','log10PVALUE']
		#if arg['--ref-genotype'] != 'miss' or 'Pr' in inf_s or 'Pd' in inf_s:
		header = ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','SNPidx1','SNPidx2','DELTA','ED','G','PVALUE','log10PVALUE']
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
			for p,gen in genotype.items():
				if '.' in gen:
					return 0,0,0
		if len(res['R']) > 2:
			fields, pools, genotype = triAllelicSites([z[0], z[1], z[3]],res,genotype)
			return fields, pools, genotype
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
		flag.append(genotype_filter(arg, genotype))
		if arg['dp_filter']:
			if dom <= min_dp or dom >= max_dp or rec <= min_dp or rec >= max_dp: #TODO
				flag.append(False)
		if arg['--het-filter']:
			#flag.append(het_filter(arg,a,b))
			flag.append(het_filter2(arg,genotype,a,b))
		if arg['--EMS']:
			flag.append(filter_EMS(arg, REF, ALT))
		if arg['--skip-indels']:
			flag.append(indels_filter(REF,ALT))
		if arg['--parental-filter'] and ('Pr' in inf_s or 'Pd' in inf_s or 'Wr' in inf_s or 'Wd' in inf_s):
			#e,f = p_al_count[0], p_al_count[1]
			#flag.append(isogenic_filter(arg,c,d,e,f))
			flag.append(parental_filter(arg,genotype))

	elif 'D' not in inf_s:
		REF,ALT,c,d = al_count[0], al_count[1], al_count[2], al_count[3]
		flag.append(genotype_filter(arg, genotype))
		if arg['dp_filter']:
			if c+d <= min_dp or c+d >= max_dp:
				flag.append(False)
		if arg['--EMS'] and (arg['--ref-genotype'] != 'miss' or 'Pr' in inf_s or 'Pd' in inf_s):
			flag.append(filter_EMS(arg, REF, ALT))
		if arg['--parental-filter'] and ('Pr' in inf_s or 'Pd' in inf_s or 'Wr' in inf_s or 'Wd' in inf_s):
			#e,f = p_al_count[0], p_al_count[1]
			#flag.append(isogenic_filter(arg,c,d,e,f))
			flag.append(parental_filter(arg,genotype))
		if arg['--skip-indels']:
			flag.append(indels_filter(REF,ALT))
	if False in flag:
		return False
	else:
		return True

def filter_qtl(arg, al_count, genotype):
	flag = [True]
	inf_s = set(arg['--data'])
	min_dp = arg['--min-depth']
	max_dp = arg['--max-depth']
	a,b,c,d = int(al_count[2]), int(al_count[3]), int(al_count[4]), int(al_count[5])
	high = a + b
	low = c + d
	dp = high + low
	GT_H = set(genotype['D'])
	GT_L = set(genotype['R'])
	flag.append(False if GT_H == GT_L == {'1'} else True)
	if arg['dp_filter']:
		flag.append(True if dp <= max_dp and dp >= min_dp else False)
	#flag.append(False if b/(a+b) == d/(c+d) else True)
	if arg['--skip-indels']:
		flag.append(indels_filter(al_count[0],al_count[1]))
	if False in flag:
		return False
	else:
		return True

def indels_filter(ref, alt):
	if len(ref) > 1 or len(alt) > 1:
		return False
	elif len(ref) == 1 and len(alt) == 1:
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

def het_filter2(arg, genotype, a, b):
	min_ratio = arg['--min-ratio']
	max_ratio = arg['--max-ratio']
	GT_dom = set(genotype['D'])
	dom = a + b
	flag = list()
	flag.append(False if GT_dom != {'0','1'} else True)
	if dom > 0:
		flag.append(True if (a/dom) <= max_ratio and (a/dom) >= min_ratio else False)
	if False in flag:
		return False
	else:
		return True


def isogenic_filter(arg,c,d,e,f):
	min_dp = arg['--min-depth']
	min_ratio = arg['--min-ratio']
	max_ratio = arg['--max-ratio']
	#min_error = arg['--min-error']
	if c + d >= min_dp and e + f >= min_dp and  min_ratio*c < min_dp < max_ratio*c and min_ratio*e < min_dp < max_ratio*e:
		return False
	#elif c == c+d and e >= min_error:
	#	return False
	#else: 
	#	return True
	
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
		if arg['--mutant-pool'] == 'R' and 'Wr' in inf_s:
			flag.append(False if GT_par == {'1'} else True)
		if 'Pr' in inf_s and arg['--mutant-pool'] == 'R':
			flag.append(False if GT_par == {'0'} else True)
		if arg['--mutant-pool'] == 'D':
			flag.append(False if '0' in GT_par else True)
			if '0' in GT_par:
				print('parental filter', arg['poss'], arg['count'],arg['pcount'], arg['genot'])
	#elif 'Wr' in inf_s or 'Wd' in inf_s:
	#	flag.append(False if GT_par == {'0','1'} or GT_par == {'1'} else True)
	elif ('Pd' in inf_s or 'Wd' in inf_s) and arg['--mutant-pool'] == 'R':
		flag.append(False if '1' in GT_par else True)
		if '1' in GT_par:
			print('parental filter', arg['poss'], arg['count'],arg['pcount'], arg['genot'])
		#TODO-flag.append(False if GT_par == GT_rec else True)
	elif ('Pd' in inf_s or 'Wd' in inf_s) and arg['--mutant-pool'] == 'D':
		if 'Wd' in inf_s:
			flag.append(False if '0' in GT_par else True)
		if 'Pd' in inf_s:
			flag.append(False if GT_par == {'1'} else False)
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

def triAllelicSites(fields, pools, genotype):
	ref = fields[2]		#alelo de la referencia lineal
	gens = list()
	for p,g in genotype.items():
		gens += g
	if '0' in gens or len(pools['R']) > 3:
		return 0,0,0
	else:
		for p,c in pools.items():
			del c[ref]
		new_ref = list(pools['R'].keys())[0]
		fields[2] = new_ref
		new_genotype = dict()
		translate = {'1':'0', '2':'1'}
		for p,g in genotype.items():
			g = [translate[i] for i in g]
			new_genotype[p] = g
		return fields, pools, new_genotype
	
def normalize(pools, REF, arg, genotype, r_min=0.03):
	data = arg['data_w']
	inf_s = set(data)
	if 'qtl' in arg.keys():
		data = arg['--data']
		inf_s = set(data)
	wt = arg['wt']
	ref = arg['--ref-genotype']
	rec = pools['R']
	alle = [a for a,v in rec.items()]
	#REF_idx = alle.index(REF)
	arg['reorder'] = 0
	al_count = 0
	wt_l = 0
	if len(alle) > 2:
		return 0,0,0
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
			if p_dom[p_al[0]] > 0 and p_dom[p_al[1]]/(p_dom[p_al[0]] + p_dom[p_al[1]]) < 0.03:# or len(p_al) > 2:
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
				al_count, wt_l, genotype = 0,0,0
	elif len(inf_s) == 3:
		if 'Pr' in inf_s:
			parental = pools['Pr']
		elif 'Pd' in inf_s:
			parental = pools['Pd']
		p_al = sorted(parental, key=lambda key: parental[key], reverse=True)
		dom = pools['D']
		r_min = 0.03
		if parental[p_al[0]] > 0 and parental[p_al[1]]/(parental[p_al[0]] + parental[p_al[1]]) < r_min:	# Asegurar que el parental esta en homocigosis
			p_al2 = [key for key,arg in parental.items()]
			if 'Pd' in inf_s:
				arg['reorder'] = 1 if p_al != p_al2 else 0
				a = dom[p_al[0]]
				b = dom[p_al[1]]
				c = rec[p_al[0]]
				d = rec[p_al[1]]
				e = parental[p_al[0]]
				f = parental[p_al[1]]
				al_count = [p_al[0], p_al[1], a, b, c, d]
				wt_l = [e, f] 
			elif 'Pr' in inf_s:
				arg['reorder'] = 0 if p_al != p_al2 else 1
				a = dom[p_al[1]]
				b = dom[p_al[0]]
				c = rec[p_al[1]]
				d = rec[p_al[0]]
				e = parental[p_al[1]]
				f = parental[p_al[0]]
				al_count = [p_al[1], p_al[0], a, b, c, d]
				wt_l = [e, f]
		else:
			al_count, wt_l, genotype = 0, 0, 0
	if wt:
		if arg['reorder'] == 1:
			wt_l = [pools[wt][alle[1]], pools[wt][alle[0]]]
		elif arg['reorder'] == 0:
			wt_l = [pools[wt][alle[0]], pools[wt][alle[1]]]
	if arg['reorder'] == 1:
		new_genotype = dict()
		translate = {'1':'0', '0':'1'}
		for p, gen in genotype.items():
			gen = [translate[i] for i in gen]
			new_genotype[p] = gen
		return al_count, wt_l, new_genotype
	new_genotype = genotype
	return al_count, wt_l, new_genotype
	

def check_annotate_args(arg):
	arg['--transl-table'] = int(arg['--transl-table'])
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
			if int(Treg[0]) > int(Treg[1]):
				print('Error: Enter a valid region (-R chromName:startPos-endPos)', file=sys.stderr)
				sys.exit()
			arg['--region'] = [Tchrom,int(Treg[0]),int(Treg[1])]
		else:
			print('Error: Enter region correctly (-R chromName:startPos-endPos)', file=sys.stderr)
			sys.exit()
	if 'R' in arg['--data'] and 'D' in arg['--data']:
		header=['#CHROM','POS','DOM','REC', 'DPdom_1','DPrec_1','DPdom_2','DPrec_2','TYPE','ID','PARENT','STRAND',\
			'CODON_change','AA_change','INFO']
	else:
		header=['#CHROM','POS','DOM','REC', 'DPdom_1','DPrec_1','TYPE','ID','PARENT','STRAND',\
			'CODON_change','AA_change','INFO']
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
						print('Warning: Phase has to be an integer in CDS {}:{}-{} and was discarded.\n'.format(seqid,start,end), file=sys.stderr)
					else:			
						gff[type_].append([seqid, type_, int(start), int(end), strand, phase, dict_att['ID'], dict_att['Parent'], dict_att['Name']])
	gff['gene_structures'] = gff['CDS'] + gff['five_prime_UTR'] + gff['three_prime_UTR']
	gff['gene'] = gff['gene'] + gff['ncRNA_gene'] + gff['pseudogene']	#all of these types are considered as genes
	gff['gene'].sort(key=lambda row:row[2]) #ordered by pos
	gff['mRNA'].sort(key=lambda row:(row[6], row[2])) #ordered by id and pos
	gff['exon'].sort(key=lambda row:(row[7], row[2])) #ordered by parent and pos
	gff['CDS'].sort(key=lambda row:(row[7], row[2])) #ordered by parent and pos
	gff['five_prime_UTR'].sort(key=lambda row:(row[7], row[2])) #ordered by parent and pos
	gff['three_prime_UTR'].sort(key=lambda row:(row[7], row[2])) #ordered by parent and pos
	gff['gene_structures'].sort(key=lambda row:(row[7], row[2]))
	ordered = ['gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR','gene_structures']
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


def find_effect(before_nt,after_nt,strand, arg):
	before_nt = Seq(before_nt)
	after_nt = Seq(after_nt)
	tab = arg['--transl-table']
	#:
	#	n = 3 - len(after_nt)%3
	#	after_N = after_nt + Seq('N')*n
	if strand == '-':
		before_nt = before_nt.reverse_complement()
		after_nt = after_nt.reverse_complement()
	if arg['variant'] == 'insertion':
		n = len(after_nt) % 3
		after = after_nt[:len(after_nt)-n].translate(table=tab)
		before = before_nt.translate(table=tab)
	elif arg['variant'] == 'deletion':
		n = len(before_nt) % 3
		before = before_nt[:len(before_nt)-n].translate(table=tab)
		after = after_nt.translate(table=tab)
	elif arg['variant'] == 'substitution':
		before = before_nt.translate(table=tab)
		after = after_nt.translate(table=tab)
	change = list()
	if (before[0] == after[0]) and len(before_nt) == len(after_nt):
		change.append('Synonymous:.')
	elif (before[0] != after[0]):
		if (after[0] == '*'):
			change.append('Nonsynonymous:nonsense')
		elif (before[0] == '*'):
			change.append('Nonsynonymous:nonstop')
		else:
			change.append('Nonsynonymous:missense')
	if len(before_nt) < len(after_nt):
		if (len(after_nt) - len(before_nt))%3 == 0:
			change.append('Insertion:in-frame')
		else:
			change.append('Insertion:frameshift')
	return [before, after, before_nt, after_nt, ';'.join(change)]

def filter_EMS(arg, REF, ALT):
	if len(REF) > 1 or len(ALT) > 1:
		return True
	if arg['--mutant-pool'] == 'R':
		if (REF == 'G' and ALT == 'A') or (REF == 'C' and ALT == 'T'):
			return True
		else:
			return False
	elif arg['--mutant-pool'] == 'D':
		if (REF == 'A' and ALT == 'G') or (REF == 'T' and ALT == 'C'):	
			return True
		else:
			#print('EMS filter',arg['poss'], arg['count'], arg['pcount'], arg['genot'])
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
	result['CODON_ref'], result['CODON_alt'], result['AA_ref'],result['AA_alt']= '.','.','.','.'
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

def is_indel(row,arg):
	reorder = row['reorder']
	ref = row['DOM'] if reorder == 0 else row['REC']
	alt = row['REC'] if reorder == 0 else row['DOM']
	if len(alt) > len(ref):
		ref,alt=NWSellers(ref, alt, 0, -1)
		row['DOM'] = ref
		row['REC'] = alt
		arg['variant'] = 'insertion'
		check_mutation2(row, arg)
	elif len(alt) < len(ref):
		
		alt,ref=NWSellers(alt, ref, 0, -1)
		row['DOM'] = ref
		row['REC'] = alt
		arg['variant'] = 'deletion'
		check_deletion(row, arg, ref, alt)
	arg['indel'] = list()


def check_deletion(row, arg, ref, alt):
	result = {h:row[h] for h in arg['header2'] if h in arg['header2'] and h in arg['header']}
	result2 = {h:'.' for h in arg['header'] if h not in arg['header2']}
	result2['INFO'] = dict()
	result.update(result2)
	pos = int(row['POS'])
	coor_i = int(pos)
	coor_f = int(pos) + len(ref) - 1
	result['CODON_ref'], result['CODON_alt'], result['AA_ref'], result['AA_alt'] = '.','.','.','.'
	gff = arg['gff']
	bi,ei = find_row(gff['gene'], coor_i, 0, len(gff['gene']))
	bf,ef = find_row(gff['gene'], coor_f, 0, len(gff['gene']))

	if ei < bi and ef < bf and ei == ef and bi == bf: #INTERGENIC, BETWEEN GENES
		if len(gff['gene']) == bi:
			dis1 = arg['--contigs'][arg['--region'][0]] - pos
			result['INFO']['right'] = 'chromosome_end'+':'+ str(dis1)
		else:
			dis1 = gff['gene'][bi][2] - coor_f
			result['INFO']['right'] = gff['gene'][bi][6] +':'+ str(dis1)
		dis2 = coor_i - gff['gene'][ei][3]
		result['TYPE'] = 'intergenic'
		result['INFO']['left'] = gff['gene'][ei][6] +':'+ str(dis2)
		write_annotate_line(arg, result, row)
		result['INFO'] = dict()
	
	elif ei < bi and ef < bf and ei != ef and bi != bf:	#INTERGENIC, COMPLETE GENES DELETION
		genes_deleted = [gene[6] for gene in gff['gene'][bi:ef+1]]
		result['TYPE'] = 'multi_gene'
		result['INFO']['effect'] = 'genes_deleted:' + ','.join(genes_deleted)
		write_annotate_line(arg, result, row)
		result['INFO'] = dict()

	elif bi > ei and bf == ef:	#LEFT INTERGENIC, RIGHT GENIC
		if len(gff['gene'][bi:ef+1]) > 1:
			genes_deleted = gff['gene'][bi:ef+1]
			c_deleted = [g[6] for g in genes_deleted[:-1]]
			gene = genes_deleted[-1] 
			result['TYPE'] = 'multi_gene'
			size_deleted = coor_f - gene[2]
			if gene[4] == '+':
				result['INFO']['5_prime_deletion'] = gene[6]+':'+str(size_deleted)
			else:
				result['INFO']['3_prime_deletion'] = gene[6]+':'+str(size_deleted)
			
			result['INFO']['effect'] = 'genes_deleted:' + ','.join(c_deleted)
			write_annotate_line(arg, result, row)
			result['INFO'] = dict()
		else:
			gene = gff['gene'][bi:ef+1][0]
			result['TYPE'] = 'gene'
			result['ID'] = gene[6]
			result['STRAND'] = gene[4]
			size_deleted = coor_f - gene[2]
			if result['STRAND'] == '+':
				result['INFO']['5_prime_deletion'] = gene[6] + ':' + str(size_deleted)
				write_annotate_line(arg, result, row)
			else:
				result['INFO']['3_prime_deletion'] = gene[6] + ':' + str(size_deleted)
				write_annotate_line(arg, result, row)
			result['INFO'] = dict()

	elif bf > ef and bi == ei: #INTERGENIC RIGHT, GENIC LEFT
		if len(gff['gene'][bi:ef+1]) > 1:
			genes_deleted = gff['gene'][bi:ef+1]
			c_deleted = [g[6] for g in genes_deleted[1:]]
			result['TYPE'] = 'multi_gene'
			gene = genes_deleted[0]
			size_deleted = gene[3] - coor_i
			if gene[4] == '+':
				result['INFO']['3_prime_deletion'] = gene[6]+':'+str(size_deleted)
			else:
				result['INFO']['5_prime_deletion'] = gene[6]+':'+str(size_deleted)

			result['INFO']['effect'] = 'genes_deleted:' + ','.join(c_deleted)
			write_annotate_line(arg, result, row)
			result['INFO'] = dict()
		else:
			gene = gff['gene'][bi:ef+1][0]
			result['TYPE'] = 'gene'
			result['ID'] = gene[6]
			result['STRAND'] = gene[4]
			size_deleted = gene[3] - coor_i
			if result['STRAND'] == '+':
				result['INFO']['3_prime_deletion'] = gene[6]+':'+ str(size_deleted)
				write_annotate_line(arg, result, row)
			else:
				result['INFO']['5_prime_deletion'] = gene[6]+':'+ str(size_deleted)
				write_annotate_line(arg, result, row)
			result['INFO'] = dict()

	elif bf == ef and bi == ei:
		genes = gff['gene'][bi:ef+1]
		if len(genes) > 1:
			genes_deleted = [gene[6] for gene in genes[1:-1]]
			result['TYPE'] = 'multi_gene'
			gene_left = genes[0]
			gene_right = genes[-1]
			size_deleted_l = gene_left[3] - coor_i
			size_deleted_r = coor_f - gene_right[2]
			if gene_left[4] == '+':
				result['INFO']['3_prime_deletion'] = gene_left[6]+':'+str(size_deleted_l)
			elif gene_left[4] == '-':
				result['INFO']['5_prime_deletion'] = gene_left[6]+':'+str(size_deleted_l)

			if gene_right[4] == '+':
				result['INFO']['5_prime_deletion'] = gene_right[6]+':'+str(size_deleted_r)
			elif gene_right[4] == '-':
				result['INFO']['3_prime_deletion'] = gene_right[6]+':'+str(size_deleted_r)
			
			if genes_deleted:
				result['INFO']['effect'] = 'genes_deleted:'+','.join(genes_deleted)
			write_annotate_line(arg, result, row)
			result['INFO'] = dict()

		else:
			gene = genes[0]
			gene_id = gene[6]
			result['TYPE'] = 'gene'
			result['STRAND'] = gene[4]
			b1,e1 = find_row_name(gff['mRNA'], gene_id, 7)
			b1i,e1i = find_row(gff['mRNA'], coor_i, b1, e1+1)
			b1f,e1f = find_row(gff['mRNA'], coor_f, b1, e1+1)
			if gff['mRNA'][b1i:e1f+1]:
				for mRNA in gff['mRNA'][b1i:e1f+1]:
					if coor_i <= mRNA[3] and coor_f >= mRNA[2]:	#check if deletion is in range for this iso-form
						genic_deletion(arg, result, row, coor_i, coor_f, mRNA)
			else:
				find_nc_gene(arg, gff, gene_id, pos, result, row)

	elif bi == bf and ei == ef and bi < ef: ##genic, in + and - strand
		genes = gff['gene'][bi:ef+1]
		for gene in genes:
			gene_id = gene[6]
			result['TYPE'] = 'gene'
			result['STRAND'] = gene[4]
			b1,e1 = find_row_name(gff['mRNA'], gene_id, 7)
			b1i,e1i = find_row(gff['mRNA'], coor_i, b1, e1+1)
			b1f,e1f = find_row(gff['mRNA'], coor_f, b1, e1+1)
			if gff['mRNA'][b1i:e1f+1]:
				for mRNA in gff['mRNA'][b1i:e1f+1]:
					if coor_i <= mRNA[3] and coor_f >= mRNA[2]:	#check if deletion is in range for this iso-form
						genic_deletion(arg, result, row, coor_i, coor_f, mRNA)
			else:
				find_nc_gene(arg, gff, gene_id, pos, result, row)

def genic_deletion(arg, result, row, coor_i, coor_f, mRNA):
	gff = arg['gff']
	mRNA_id = mRNA[6]
	b,e = find_row_name(gff['exon'], mRNA_id, 7)
	bi,ei = find_row(gff['exon'], coor_i, b, e+1)
	bf,ef = find_row(gff['exon'], coor_f, b, e+1)
	result['PARENT'] = mRNA_id

	if ei < bi and ef < bf and ei == ef and bi == bf:	#Complete intronic deletion
		result['TYPE'] = 'intron'
		dis1 = gff['exon'][bi][2] - coor_f
		dis2 = coor_i - gff['exon'][ei][3]
		result['INFO']['left'] = gff['exon'][ei][8] +':'+ str(dis2)
		result['INFO']['right'] = gff['exon'][bi][8] +':'+ str(dis1)
		if gff['exon'][ei][7] == gff['exon'][bi][7]:
			if dis2 in {1,2,3}:
				if gff['exon'][ei][4] == '+':
					result['INFO']['5_splice_site'] = 'intron-boundary:'+gff['exon'][ei][8]
				else:
					result['INFO']['3_splice_site'] = 'intron-boundary:'+gff['exon'][ei][8]
			if dis1 in {1,2,3}:
				if gff['exon'][bi][4] == '+':
					result['INFO']['3_splice_site'] = 'intron-boundary:'+gff['exon'][bi][8]
				else:
					result['INFO']['5_splice_site'] = 'intron-boundary:'+gff['exon'][bi][8]
			write_annotate_line(arg, result, row)
			result['INFO'] = dict()
	elif ei < bi and ef < bf and ei != ef and bi != bf:	#complete exon deletion
		exons_deleted = [exon[8] for exon in gff['exon'][bi:ef+1]]
		result['TYPE'] = 'multi_exon'
		result['INFO']['effect'] = 'exons_deleted:' + ','.join(exons_deleted)
		write_annotate_line(arg, result, row)
		result['INFO'] = dict()

	elif bi > ei and bf == ef:	#LEFT INTRON, RIGHT EXON
		if len(gff['exon'][bi:ef+1]) > 1:
			exons_deleted = gff['exon'][bi:ef+1]
			c_deleted = [e[8] for e in exons_deleted[:-1]]
			exon = exons_deleted[-1]
			result['TYPE'] = 'multi_exon'
			result['STRAND'] = exon[4]
			if gff['exon'][bi-1][7] == exon[7]:
				if exon[4] == '+':
					result['INFO']['3_splice_site_deletion'] = 'intron/exon-boundary:'+exon[8]
				else:
					result['INFO']['5_splice_site_deletion'] = 'exon/intron-boundary:'+exon[8]

			result['INFO']['effect'] = 'exons_deleted:' + ','.join(c_deleted)
			write_annotate_line(arg, result, row)
			result['INFO'] = dict()
		elif len(gff['exon'][bi:ef+1]) == 1:
			exon = gff['exon'][bi:ef+1][0]
			result['TYPE'] = 'exon'
			result['STRAND'] = exon[4]
			if gff['exon'][bi-1][7] == exon[7]:
				if exon[4] == '+':
					result['INFO']['3_splice_site_deletion'] = 'intron/exon-boundary:'+exon[8]
				else:
					result['INFO']['5_splice_site_deletion'] = 'exon/intron-boundary:'+exon[8]

				write_annotate_line(arg, result, row)
				result['INFO'] = dict()

	elif bf > ef and bi == ei:	#LEFT EXON, RIGHT INTRON
		if len(gff['exon'][bi:ef+1]) > 1:
			exons_deleted = gff['exon'][bi:ef+1]
			c_deleted = [e[8] for e in exons_deleted[1:]]
			exon = exons_deleted[0]
			result['TYPE'] = 'multi_exon'
			result['STRAND'] = exon[4]
			if gff['exon'][bi+1][7] == exon[7]:
				if exon[4] == '+':
					result['INFO']['5_splice_site_deletion'] = 'exon/intron-boundary:'+exon[8]
				else:
					result['INFO']['3_splice_site_deletion'] = 'intron/exon-boundary:'+exon[8]

			result['INFO']['effect'] = 'exons_deleted:' + ','.join(c_deleted)
			write_annotate_line(arg, result, row)
			result['INFO'] = dict()
		elif len(gff['exon'][bi:ef+1]) == 1:
			exon = gff['exon'][bi:ef+1][0]
			result['TYPE'] = 'exon'
			result['STRAND'] = exon[4]
			if gff['exon'][bi+1][7] == exon[7]:
				if exon[4] == '+':
					result['INFO']['5_splice_site_deletion'] = 'exon/intron-boundary:'+exon[8]
				else:
					result['INFO']['3_splice_site_deletion'] = 'intron/exon-boundary:'+exon[8]

				write_annotate_line(arg, result, row)
				result['INFO'] = dict()
	elif bf == ef and bi == ei:
		if bi != bf:
			exons = gff['exon'][bi:ef+1]
			exons_deleted = [exon[8] for exon in exons[1:-1]]
			result['TYPE'] = 'multi_exon'
			exon_left = exons[0]
			exon_right = exons[-1]
			#size_deleted_l = gene_left[3] - coor_i
			#size_deleted_r = coor_f - gene_right[2]
			if gff['exon'][bi+1][7] == exon_left[7]:
				if exon_left[4] == '+':
					result['INFO']['5_splice_site_deletion'] ='exon/intron-boundary:'+ exon_left[8]
				elif exon_left[4] == '-':
					result['INFO']['3_splice_site_deletion'] ='intron/exon-boundary:'+ exon_left[8]
			if gff['exon'][ef-1][7] == exon_right[7]:
				if exon_right[4] == '+':
					result['INFO']['3_splice_site_deletion'] ='intron/exon-boundary:'+ exon_right[8]
				elif exon_right[4] == '-':
					result['INFO']['5_splice_site_deletion'] ='exon/intron-boundary:'+ exon_right[8]

			if exons_deleted:
				result['INFO']['effect'] = 'exons_deleted:'+','.join(exons_deleted)
			write_annotate_line(arg, result, row)
			result['INFO'] = dict()
		elif bi == bf:
			b2, e2 = find_row_name(gff['gene_structures'], mRNA_id, 7)
			b2i, e2i = find_row(gff['gene_structures'], coor_i, b2, e2+1)
			b2f, e2f = find_row(gff['gene_structures'], coor_f, b2, e2+1)
			if gff['gene_structures'][b2i:e2i+1] and gff['gene_structures'][b2f:e2f+1]:
				left = gff['gene_structures'][b2i:e2i+1][0]
				right = gff['gene_structures'][b2f:e2f+1][0]
				if left == right and left[1] == 'three_prime_UTR':
					result['TYPE'] = 'three_prime_UTR'
					result['PARENT'] = mRNA_id
					write_annotate_line(arg, result, row)
					result['INFO'] = dict()
				elif left == right and left[1] == 'five_prime_UTR':
					result['TYPE'] = 'five_prime_UTR'		
					result['INFO']['effect'] = '.:.:.'
					result['PARENT'] = mRNA_id
					write_annotate_line(arg, result, row)
					result['INFO'] = dict()
				elif left == right and left[1] == 'CDS':
					result['TYPE'] = 'CDS'
					deletion = coor_f - coor_i
					if deletion % 3 == 0:
						result['INFO']['effect'] ='Deletion:in-frame'
					else:
						result['INFO']['effect'] ='Deletion:frameshift'
					write_annotate_line(arg, result, row)
					result['INFO'] = dict()
				elif left[1] == 'CDS' and right[1] == 'three_prime_UTR':
					result['TYPE'] = 'gene'
					result['STRAND'] = '+'
					deletion = right[2] - coor_i #left[3] - coor_i
					if deletion % 3 == 0:
						result['INFO']['effect'] ='Deletion:in-frame'
					else:
						result['INFO']['effect'] ='Deletion:frameshift'
					result['INFO']['left'] = 'CDS'
					result['INFO']['right'] = 'three_prime_UTR'
					write_annotate_line(arg, result, row)
					result['INFO'] = dict()
				elif right[1] == 'CDS' and left[1] == 'three_prime_UTR':
					result['TYPE'] = 'gene'
					result['STRAND'] = '-'
					deletion = coor_f - left[3]# coor_f - right[2]
					if deletion % 3 == 0:
						result['INFO']['effect'] ='Deletion:in-frame'
					else:
						result['INFO']['effect'] ='Deletion:frameshift'
					result['INFO']['left'] = 'three_prime_UTR'
					result['INFO']['right'] = 'CDS'
					write_annotate_line(arg, result, row)
					result['INFO'] = dict()
				elif right[1] == 'CDS' and left[1] == 'five_prime_UTR':
					result['TYPE'] = 'gene'
					result['STRAND'] = '+'
					deletion = coor_f - left[3] #coor_f - right[2]
					if deletion % 3 == 0:
						result['INFO']['effect'] ='Deletion:in-frame'
					else:
						result['INFO']['effect'] ='Deletion:frameshift'
					result['INFO']['left'] = 'five_prime_UTR'
					result['INFO']['right'] = 'CDS'
					write_annotate_line(arg, result, row)
					result['INFO'] = dict()
				elif left[1] == 'CDS' and right[1] == 'five_prime_UTR':
					result['TYPE'] = 'gene'
					result['STRAND'] = '-'
					deletion = right[2] - coor_i #left[3] - coor_i
					if deletion % 3 == 0:
						result['INFO']['effect'] ='Deletion:in-frame'
					else:
						result['INFO']['effect'] ='Deletion:frameshift'
					result['INFO']['left'] = 'CDS'
					result['INFO']['right'] = 'five_prime_UTR'
					write_annotate_line(arg, result, row)
					result['INFO'] = dict()

def find_nc_gene(arg, gff, gene_id, pos, result, row):
	nc_genes = ['miRNA', 'tRNA', 'ncRNA', 'ncRNA_gene', 'lnc_RNA', 'snoRNA', 'snRNA',\
		 'rRNA', 'pseudogene', 'pseudogenic_transcript', 'pre_miRNA', 'SRP_RNA']
	for type_ in nc_genes:
		b, e = find_row_name(gff[type_], gene_id, 7)
		b2,e2 = find_row(gff[type_], pos, b, e+1)
		check_nc_gene(arg,gff,type_, b2, e2, pos,result,row)


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
							exon = gff['exon'][b4:e4+1][0]
							if gff['exon'][b4-1][7] == exon[7]: #if prev exon exist
								if pos - exon[2] in {0,1,2}:
									if exon[4] == '+':
										result['INFO']['3_splice_site'] = 'exon-boundary:'+gff['exon'][e4][8]
									else:
										result['INFO']['5_splice_site'] = 'exon-boundary:'+gff['exon'][e4][8]

							if gff['exon'][b4+1][7] == exon[7]: #if next cds exist
								if exon[3] - pos in {0,1,2}:
									if exon[4] == '+':
										result['INFO']['5_splice_site'] = 'exon-boundary:'+gff['exon'][e4][8]
									else:
										result['INFO']['3_splice_site'] = 'exon-boundary:'+gff['exon'][e4][8]
				
							b5,e5 = find_row_name(gff['CDS'], mRNA_id, 7)
							b6,e6 = find_row(gff['CDS'], pos, b5, e5+1)
							b7,e7 = find_row_name(gff['five_prime_UTR'], mRNA_id, 7)
							b8,e8 = find_row(gff['five_prime_UTR'], pos, b7, e7+1)
							b9,e9 = find_row_name(gff['three_prime_UTR'], mRNA_id, 7)
							b10,e10 = find_row(gff['three_prime_UTR'], pos, b9, e9+1)
							if gff['CDS'][b6:e6+1]:
								for cds in gff['CDS'][b6:e6+1]:
									cds_id = cds[6]
									result['TYPE'] = 'CDS'
									beg, end = codon_coords(pos,cds[2],cds[3],cds[4],int(cds[5])) #target pos, cds start, cds end, cds strand, cds phase

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
										aa_before, aa_after, codon_bef, codon_aft, change = find_effect(codon_ref, codon_alt, cds[4], arg)
										result['CODON_ref'],result['CODON_alt'],result['AA_ref'],result['AA_alt'] = codon_bef, codon_aft, aa_before, aa_after
										result['INFO']['effect'] = change
										result['PARENT'] = mRNA_id
										write_annotate_line(arg, result, row)
										result['INFO'] = dict()
										result['CODON_ref'], result['CODON_alt'], result['AA_ref'],result['AA_alt']= '.','.','.','.'
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
										aa_before, aa_after, codon_bef, codon_aft, change = find_effect(codon_ref, codon_alt, cds[4], arg)
										result['CODON_ref'],result['CODON_alt'],result['AA_ref'],result['AA_alt'] = codon_bef, codon_aft, aa_before, aa_after
										result['PARENT'] = mRNA_id
										result['INFO']['effect'] = change
										write_annotate_line(arg, result, row)
										result['INFO'] = dict()
										result['CODON_ref'], result['CODON_alt'], result['AA_ref'],result['AA_alt']= '.','.','.','.'
									else:
										codon = arg['ref'].seq[beg-1:end]
										codon_ref = codon[:pos-beg]+row['DOM']+codon[pos-beg+1:]
										codon_alt = codon[:pos-beg]+row['REC']+codon[pos-beg+1:]
										aa_before, aa_after, codon_bef, codon_aft, change = find_effect(codon_ref, codon_alt, cds[4], arg)
										
										result['PARENT'] = mRNA_id
										result['CODON_ref'],result['CODON_alt'],result['AA_ref'],result['AA_alt'] = codon_bef, codon_aft, aa_before, aa_after
										result['INFO']['effect'] = change
										write_annotate_line(arg, result, row)
										result['INFO'] = dict()
										result['CODON_ref'], result['CODON_alt'], result['AA_ref'],result['AA_alt']= '.','.','.','.'
							elif gff['five_prime_UTR'][b8:e8+1]:
								#five =  gff['five_prime_UTR'][b8:e8+1][0]
								for five in  gff['five_prime_UTR'][b8:e8+1]:
									result['TYPE'] = 'five_prime_UTR'		
									if arg['variant'] == 'insertion':
										result['INFO']['effect'] = '.:.:.'
									elif arg['variant'] == 'substitution':
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
							if gff['exon'][e4][7] == gff['exon'][b4][7]:
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
			
			find_nc_gene(arg, gff, gene_id, pos, result, row)

	else:
		if len(gff['gene']) == b:
			dis2 = arg['--contigs'][arg['--region'][0]] - pos
			result['INFO']['right'] = 'chromosome_end'+':'+ str(dis2)
		else:
			dis1 = gff['gene'][b][2] - pos
			result['INFO']['right'] = gff['gene'][b][6] +':'+ str(dis1)
		dis2 = pos - gff['gene'][e][3]
		result['TYPE'] = 'intergenic'
		result['INFO']['left'] = gff['gene'][e][6] +':'+ str(dis2) if e != -1 else '.'
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
	di = ' > ' if arg['--mutant-pool'] == 'R' else ' < '
	if codonREF == '.' and codonALT == '.' and AAREF == '.' and AAALT == '.':
		result['CODON_change'] = '.'
		result['AA_change'] = '.'
	else:
		if row['reorder'] == 1:
			result['CODON_change'] = codonALT+di+codonREF
			result['AA_change'] = AAALT+di+AAREF
		elif row['reorder'] == 0:
			result['CODON_change'] = codonREF+di+codonALT
			result['AA_change'] = AAREF +di+AAALT
	result['INFO']['variant_type'] = arg['variant']
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
	if line.startswith('##GATKCommandLine='):
		write_line(line,fsal)
	if line.startswith('##source=HaplotypeCaller'):
		write_line(line,fsal)
	if line.startswith('##bcftools_callCommand='):
		write_line(line.split(';')[0]+'\n',fsal)
	if line.startswith('##reference='):
		if 'mbs' in arg.keys() or 'qtl' in arg.keys():
			write_line(line,fsal)       
	if line.startswith('##contig=<ID='):
		c = line.split('##contig=<ID=')[1].split(',')[0]
		s = line[line.find('<')+1:line.rfind('>')].split(',')
		id = s[0].split('=')[1]
		length = s[1].split('=')[1]
		arg['--contigs'][id] = int(length)
		arg['chromosomes'].append(c)
		write_line(line,fsal)
	if line.startswith('#CHROM'):
		bam_list = line.rstrip().split('\t')[9:]
		if len(arg['--data']) > len(bam_list):
			print('Error: The \"--data\" list does not match with the number of pools used.', file=sys.stderr)
			sys.exit()
		elif len(arg['--data']) <= len(bam_list):
			if 'qtl' in arg.keys():
				translate = {'D':'H','R':'L','Pd':'P'}
				d = [translate[i] for i in arg['--data']]
				arg['pools'] = {d[i]:bam_list[i] for i in range(len(d))}
				pools = '##bamFiles:genotypes='+','.join(':'.join((key,val)) for (key, val) in arg['pools'].items()) + '\n'
				write_line(pools, fsal)
			else:
				arg['pools'] = {arg['--data'][i]:bam_list[i] for i in range(len(arg['--data']))}
				pools = '##bamFiles:genotypes='+','.join(':'.join((key,val)) for (key, val) in arg['pools'].items()) + '\n'
				write_line(pools, fsal)
		if 'annotate' in arg.keys():
			chrom, pos_i, pos_f = arg['--region'][0], arg['--region'][1], arg['--region'][2]
			if chrom not in arg['--contigs'].keys():
				print('Error: Enter a valid chromosome name (-R chromName:startPos-endPos)', file=sys.stderr)
				sys.exit()
			maxim = min(pos_f, arg['--contigs'][chrom])
			minim = max(1, pos_i)
			arg['--region'] = [chrom, minim, maxim]


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
			if field == 'DOM' or field == 'REC' or field == 'ALT' or field == 'REF':
				continue
		if field in variable_descriptions.keys():
			write_line(variable_descriptions[field], fsal)
	if argv[0] == 'merge':
		line ='##maptools_mergeCommandINFO=\"All columns now represent the grouped values for the provided window size.\"\n'
		write_line(line, fsal)

def read_header_merge(arg:dict):
	lines = list()
	arg['--contigs'] = dict()
	with open(arg['--input'], 'r') as handle:
		for line in handle:
			if line.startswith('##maptools_'):
				lines.append(line)
			if line.startswith('##maptools_qtlCommand='):
				qtl_argv = line.split('=')[1].split(' ')
				result = list()
				for i in range(len(qtl_argv)):
					ar = qtl_argv[i]
					if ar == '-r' or ar == '--ref-genotype':
						ar_index = qtl_argv.index(ar)
						ref = qtl_argv[ar_index+1]
						result.append(True if ref != 'miss' else False)
					elif ar == '-d' or ar == '--data':
						ar_index = qtl_argv.index(ar)
						data = qtl_argv[ar_index+1].split(',')
						result.append(True if 'P' in data else False)
				#arg['--ref-genotype'] = True if True in result else False
				if not True in result:
					print('Error: A reference or a parental resequencing is needed in order to execute the \"merge\" command.', file=sys.stderr)
					sys.exit()
			if line.startswith('##maptools_mbsCommand='):
				mbs_argv = line.split('=')[1].split(' ')
				result = list()
				for i in range(len(mbs_argv)):
					ar = mbs_argv[i]
					if ar == '-r' or ar == '--ref-genotype':
						ar_index = mbs_argv.index(ar)
						ref = mbs_argv[ar_index+1]
						result.append(True if ref != 'miss' else False)
					elif ar == '-d' or ar == '--data':
						ar_index = mbs_argv.index(ar)
						data = mbs_argv[ar_index+1].split(',')
						result.append(True if 'Pd' in data or 'Pr' in data else False)
				#arg['--ref-genotype'] = True if True in result else False
				if not True in result:
					print('Error: A reference or a parental resequencing is needed in order to execute the \"merge\" command.', file=sys.stderr)
					sys.exit()
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
				s = line[line.find('<')+1:line.rfind('>')].split(',')
				id = s[0].split('=')[1]
				length = s[1].split('=')[1]
				arg['--contigs'][id] = int(length)
				lines.append(line)
			if line.startswith('##bamFiles'):
				lines.append(line)
			if line.startswith('#CHROM'):
				header = line.split('\t')
				header = [field.rstrip() for field in header]
				arg['header'] = header
				arg['header2'] = header
				if 'REF' in header and 'ALT' in header:
					translate = {'REF':'DOM','ALT':'REC','DPref_1':'DPdom_1','DPalt_1':'DPrec_1','DPref_2':'DPdom_2','DPalt_2':'DPrec_2'}
					header2 = [i if i not in translate.keys() else translate[i] for i in arg['header']]
					arg['header2'] = header2
				break
	return lines

def check_chroms(arg):
	chromosomes = arg['--chromosomes']
	for chrom in chromosomes:
		if chrom not in arg['--contigs'].keys():
			print('Error: {} is not in the contig ID list.'.format(chrom), file=sys.stderr)
			sys.exit()


#NW-Sellers
def inicializacion_NWS(n_rows,n_cols,gap):
  matriz = np.full([n_rows, n_cols], 0)
  for i in range(1, n_rows):
    matriz[i,0] = matriz[i-1, 0] + -100
  for j in range(1, n_cols):
    matriz[0,j] = matriz[0, j-1] + -100
  return matriz

def rellenado_NWS(matriz, A, B, n_rows, n_cols, gap, mismatch):
  for i in range(1, n_rows):
    for j in range(1, n_cols):
      izquierda = matriz[i,j-1] + gap
      arriba = matriz[i-1, j] + gap
      if A[j-1] == B[i-1]:
        diagonal = matriz[i-1, j-1] + 1
      else:
        diagonal = matriz[i-1, j-1] + mismatch

      matriz[i,j] = max([arriba, izquierda, diagonal])

  return matriz

def vuelta_atras_NWS(matriz, i, j, A, B, gap):
  alin_A = str()
  alin_B = str()

  while i != 0 or j != 0:

    if matriz[i,j] == matriz[i-1,j-1] + 1 and A[j-1] == B[i-1]: #desplazamiento en diagonal
      alin_A = A[j-1] + alin_A
      alin_B = B[i-1] + alin_B
      i = i - 1
      j = j - 1
    elif matriz[i,j] == matriz[i,j-1] + gap:      #desplazamiento en horizontal
      alin_A = A[j-1] + alin_A
      alin_B = "-" + alin_B
      j = j - 1
    elif matriz[i,j] == matriz[i-1, j] + gap:   #desplazamiento en vertical
      alin_A = "-" + alin_A
      alin_B = B[i-1] + alin_B
      i = i - 1


  return alin_A, alin_B

def NWSellers(A, B, gap, mismatch):
  n_cols = len(A) + 1
  n_rows = len(B) + 1
  m = inicializacion_NWS(n_rows, n_cols, gap)
  m = rellenado_NWS(m,A,B,n_rows, n_cols, gap, mismatch)
  i = n_rows - 1
  j = n_cols - 1
  a, b = vuelta_atras_NWS(m,i,j,A, B, gap)
  ref = a[0]
  alt = b[0:a.count('-')+1]
  return ref,alt	