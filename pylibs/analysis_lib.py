from logging import exception
import re
import sys
import os
from math import log
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import hypergeom
from scipy.spatial import distance
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
	g = 2*(a*log(a/n1)+b*log(b/n2)+c*log(c/n3)+d*log(d/n4))
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
		for i in range(len(data)):
			inf[data[i]] = c
			c += 1
		inf_s = set(data)
		if len(inf_s) == 1:
			data1 = z[9].split(':')
			pool = dict()
			#GT1 = data1[GT_index].split('/')
			AD = data1[AD_index].split(',')
			ref = z[3]
			alt = z[4].split(',')
			pool[ref] = int(AD[0])
			for i in range(len(alt)):
				pool[alt[i]] = int(AD[i+1])
			res[data[0]] = pool
			# #contains all the genotypes of the line
		elif len(inf_s) >= 2:
			if inf_s == {'D','R'} or inf_s == {'D','R','Pr'} or inf_s == {'D','R','Pd'}:
				data1 = z[inf['D']].split(':')
				data2 = z[inf['R']].split(':')
				GT1 = data1[GT_index].split('/')
				GT2 = data2[GT_index].split('/')
				gt = set(GT1+GT2)
				if len(gt) == 1 or '.' in gt: #if all the genotypes are equal the line is discarded
					return 0,0
				for p,c in inf.items():
					ref = z[3]
					alt = z[4].split(',')
					AD = z[c].split(':')[AD_index].split(',')
					pool = {ref:int(AD[0])}
					for i in range(len(alt)):
						pool[alt[i]] = int(AD[i+1])
					res[p] = pool
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
