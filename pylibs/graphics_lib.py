from pylibs.analysis_lib import *
from pylibs.constants import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
import scipy.stats as st
import sys
import os
from math import log
import json
from docopt import docopt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter,
                               AutoMinorLocator)

# warnings.filterwarnings("ignore")
# sys.path.append(os.getcwd())
pd.options.mode.chained_assignment = None
# ['#CHROM','POS','REF','REC','DPdom_1','DPrec_1','MAX_SNPidx2','BOOST']
# ['#CHROM','POS','REF','REC','DPdom_1','DPrec_1','SNPidx1','BOOST']
# ['#CHROM','POS','REF','REC','DPdom_1','DPrec_1','DPdom_2','DPrec_2','SNPidx1','SNPidx2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE']
# ['#CHROM','POS','REF','REC','DPdom_1','DPrec_1','DPdom_2','DPrec_2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE']
# ['#CHROM','POS','REF','REC','DPdom_1','DPrec_1','DPdom_2','DPrec_2','PVALUE','log10PVALUE']
# ['#CHROM','POS','REF','REC','DPdom_1','DPrec_1','DPdom_2','DPrec_2','SNPidx1','SNPidx2','DELTA','PVALUE','log10PVALUE']

def read_header_plot(arg):
    arg['--contigs'] = dict()
    with open(arg['--input'],'r') as handle:
        for line in handle:
            if line.startswith('##maptools_mbsCommand='):
                if 'mbsplot' not in arg.keys():
                    print('Error: the input data must come from the qtl command', file=sys.stderr)
                    sys.exit()
            if line.startswith('##maptools_mergeCommand='):
                merge_argv = line.split('=')[1].split(' ')
                if '-w' or '--window' in merge_argv:
                    names = ['-w','--window']
                    for w in names:
                        if w in merge_argv:
                            w_index = merge_argv.index(w)
                            arg['--window'] = merge_argv[w_index+1]
                if '-D' or '--distance' in merge_argv:
                    names = ['-D','--distance']
                    for w in names:
                        if w in merge_argv:
                            w_index = merge_argv.index(w)
                            arg['--distance'] = merge_argv[w_index+1]
            if line.startswith('##maptools_qtlCommand='):
                if 'qtlplot' not in arg.keys():
                    print('Error: the input data must come from the mbs command', file=sys.stderr)
                    sys.exit()
                names = ['-r', '--ref-genotype','-d','--data']
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

                arg['--ref-genotype'] = True if True in result else False

            if line.startswith('##contig=<'):
                s = line[line.find('<')+1:line.rfind('>')].split(',')
                id = s[0].split('=')[1]
                length = s[1].split('=')[1]
                arg['--contigs'][id] = int(length)

            if line.startswith('#CHROM'):
                arg['header'] = line.rstrip().split('\t')
                if 'qtlplot' in arg.keys():
                    translate = {'REF':'DOM','ALT':'REC','DPref_1':'DPdom_1','DPalt_1':'DPrec_1','DPref_2':'DPdom_2','DPalt_2':'DPrec_2'}
                    header2 = [i if i not in translate.keys() else translate[i] for i in arg['header']]
                    arg['header'] = header2
                break
    return arg

                
def test_plot(arg, __doc__):
    families = sorted(set([f.name for f in fm.fontManager.ttflist]))
    fonts = ["Arial","Liberation Sans", "Ubuntu", "DejaVu Sans", "Ubuntu", "Free Sans", "Droid Sans Fallback"]
    for font in fonts:
        if font in families:
            break
    plt.rcParams['font.family'] = font
    global fields
    inp_f = arg['--input']
    if not inp_f:
        print(__doc__, end='\n', file=sys.stderr)
        sys.exit()
    try:
        f = open(inp_f, 'r')
    except FileNotFoundError:
        print('Error: The input file {} does not exist'.format(inp_f), file=sys.stderr)
        sys.exit()
    wd = os.getcwd()
    try:
        os.makedirs(wd+'/'+arg['--outdir'])
    except FileExistsError:
        print('Warning: the output directory already exists', file=sys.stderr)
        pass
    arg['--outdir'] = wd+'/'+arg['--outdir']+'/'
    if arg['--captions']:
        arg['captions_dir'] = arg['--outdir']+'captions/'
        try:
            os.makedirs(arg['captions_dir'])
        except FileExistsError:
            pass

    read_palette(arg)

    if arg['--multi-chrom'] == True and len(arg['--chromosomes']) == 1:
        arg['--multi-chrom'] = False
    arg['--alpha'] = float(arg['--alpha'])
    arg['lim'] = 1e-90
    arg['--output-type'] = '.'+arg['--output-type']

    if arg['--moving-avg'] == False or arg['--moving-avg'] == None :
        arg['--moving-avg'] = False
    else:
        arg['--moving-avg'] = int(arg['--moving-avg'])
        if arg['--moving-avg'] <= 0:
            print(
                'Error: The window size for moving average (-A) must be an integer higher than zero.', file=sys.stderr)
            sys.exit()
    if arg['--distance-avg'] == False or arg['--distance-avg'] == None:
        arg['--distance-avg'] = False
    else:
        arg['--distance-avg'] = int(arg['--distance-avg'])
        if arg['--distance-avg'] <= 0:
            print(
                'Error: The window size for distance average (-W) must be an integer higher than zero.', file=sys.stderr)
            sys.exit()
    if arg['--moving-avg'] != False and arg['--distance-avg'] != False:
        #print(arg['--moving-avg'], arg['--distance-avg'])
        print('Error: Options -A and -W are incompatible', file=sys.stderr)
        sys.exit()
    return arg

def load_dataframe_plotting(arg):
    inp_f = arg['--input']
    inp_ext = inp_f.split('.')[-1]
    if inp_ext == 'csv':
        sep_ = ','
    else:
        sep_ = '\t'
    df = pd.read_csv(inp_f, sep=sep_, dtype=fields, names=arg['header'], comment='#')
    chroms = df['#CHROM'].unique()
    if arg['--chromosomes'] == 'all':
        arg['--chromosomes'] = list(chroms)
    else:
        chroms_ = arg['--chromosomes'].split(',')
        arg['--chromosomes'] = chroms_
    check_chroms(arg)
    arg['labs'] = {arg['--chromosomes'][i]:chr(97+i) for i in range(len(arg['--chromosomes']))}
    chrom_lists = [arg['--chromosomes'][i:i+6] for i in range(0, len(arg['--chromosomes']), 6)]
    arg['chrom_lists'] = chrom_lists
    arg['--fields'] = list(df.columns)
    arg['contigs'] = {ch:arg['--contigs'][ch] for ch in arg['--chromosomes']}
    arg['max_lenght'] = max(arg['contigs'].values())
    if 'log10PVALUE' in arg['--fields']:
        r = pd.DataFrame(dtype=np.float32)
        r['log10PVALUE'] = df['log10PVALUE']
        r = r.replace({-np.inf: r[np.isfinite(r)].min().min()})
        df['log10PVALUE'] = r['log10PVALUE']
    arg['n_markers'] = len(df)
        # Checking graphic types
    if 'mbsplot' in arg.keys():
        arg = check_mbs_opts(arg)
    if 'qtlplot' in arg.keys():
        arg = check_qtl_opts(arg, df)
    return arg, df

def load_dataframe(arg):
    inp_f = arg['--input']
    inp_ext = inp_f.split('.')[-1]
    if inp_ext == 'csv':
        sep_ = ','
    else:
        sep_ = '\t'
    df = pd.read_csv(inp_f, sep=sep_, dtype=fields, names=arg['header2'], comment='#')
    arg['--fields'] = list(df.columns)
    if 'DOM' in arg['--fields'] and 'REC' in arg['--fields']:
        arg['--fields'].remove('DOM')
        arg['--fields'].remove('REC')
    #elif 'REF' in arg['--fields'] and 'ALT' in arg['--fields']:
    #    arg['--fields'].remove('REF')
    #    arg['--fields'].remove('ALT')
    arg['header']=arg['--fields']
    if arg['--chromosomes'] == 'all':
        arg['--chromosomes'] = list(df['#CHROM'].unique())
    else:
        arg['--chromosomes'] = arg['--chromosomes'].split(',')
    check_chroms(arg)
    return arg, df
    
def check_merge(arg,__doc__):
    global fields
    inp_f = arg['--input']
    if not inp_f:
        print(__doc__,end='\n', file=sys.stderr)
        sys.exit()

    try:
        f = open(arg['--input'], 'r')
    except FileNotFoundError:
        print('Error: The input file {} does not exist'.format(inp_f), file=sys.stderr)
        sys.exit()
    if arg['--output-type'] in {'csv','txt'}:
        if arg['--output-type'] == 'csv':
            arg['spacer'] = ','
        elif arg['--output-type'] == 'txt':
            arg['spacer'] = '\t'
        arg['--output-type'] = '.' + arg['--output-type']
    else:
        print('Error: select a valid format.', file=sys.stderr)
        sys.exit()
    if arg['--output'] != None:
        wd = os.getcwd()
        outdir_list = arg['--output'].split('/')[:-1]
        outdir = '/'.join(outdir_list)+'/'
        try:
            os.makedirs(wd+'/'+outdir)
        except FileExistsError:
            print('Warning: the output directory already exists', file=sys.stderr)
            pass
        arg['filename'] = arg['--output'].split('/')[-1]
        arg['outdir'] = wd+'/'+outdir
    
        arg['--output'] = arg['outdir']+arg['filename']
        if arg['--output'] != None:
            arg['--output'] = check_save_an(arg, arg['filename'])
        else:
            arg['--output'] = None

    arg['lim'] = 1e-90

    if arg['--window'] != None and arg['--distance'] != None:
        print('Error: Optios -w and -D are incompatible.', file=sys.stderr)
        sys.exit()
    if arg['--window'] == None:
        arg['--window'] = False
    else:
        arg['--window'] = int(arg['--window'])
        if arg['--window'] <= 0:
            print('Error: The window size must be an integer greater than zero.', file=sys.stderr)
            sys.exit()
    
    if arg['--distance'] == None:
        arg['--distance'] = False
    else:
        arg['--distance'] = int(arg['--distance'])
        if arg['--distance'] <= 0:
            print('Error: The distance size must be an integer greater than zero.', file=sys.stderr)
            sys.exit()
    # Checking graphic types
    return arg

def read_palette(arg):
    try:
        with open("palette.json") as json_file:
            data = json.load(json_file)
    except FileNotFoundError:
        print('Error: Please put the file \"palette.json\" in the MAPtools folder.', file=sys.stderr)
        sys.exit()
    arg['DPI'] = data['DPI']
    ant = 'mbs' if 'mbsplot' in arg.keys() else 'qtl'
    palette = arg['--palette']
    if palette not in data[ant].keys():
        print('Error: Select a correct palette\'s name', file=sys.stderr)
        sys.exit()
    
    for key,val in data[ant][palette].items():
        if val == list():
            palette = 'standard'
            break
    arg['--palette'] = {k: v[0] for k, v in data[ant][palette].items()}
    arg['color_names'] = {k: v[1] for k, v in data[ant][palette].items()}
        

def check_mbs_opts(arg):
    if arg['--boost'] != None and arg['--distance-boost'] != None:
        print('Error: Options -b and -B are incompatible', file=sys.stderr)
        sys.exit()
    
    if arg['--boost'] == None:
        arg['--boost'] = False
    else:
        arg['--boost'] = int(arg['--boost'])
        if arg['--boost'] <= 0:
            print('Error: The window size for boost (-b) must be an integer higher than zero.', file=sys.stderr)
            sys.exit()
    
    if arg['--distance-boost'] == None:
        arg['--distance-boost'] = False
    else:
        arg['--distance-boost'] = int(arg['--distance-boost'])
        if arg['--distance-boost'] <= 0:
            print('Error: The window size for boost (-B) must be an integer higher than zero.', file=sys.stderr)
            sys.exit()
    arg['titles'] = titles_mbs
    arg['lines'] = lines_mbs

    if arg['--all']:
        if 'log10PVALUE' in arg['--fields']:
            arg['--pvalue'] = True
        else:
            arg['--pvalue'] = False
        if 'SNPidx1' in arg['--fields']:
            arg['--allele-freq-1'] = True
        if 'SNPidx2' in arg['--fields']:
            arg['--allele-freq-2'] = True
        if 'SNPidx1' in arg['--fields'] and 'SNPidx2' in arg['--fields']:
            arg['--combine'] = True
        if len(arg['--chromosomes']) > 1:
            arg['--multi-chrom'] = True
        if 'MAX_SNPidx2' in arg['--fields']:
            arg['--max-allele-freq2'] = True
        else:
            arg['--max-allele-freq2'] = False
    else:
        if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields'] and arg['--combine'] == True:
            print('Error: is not possible make combined graphics with your input data', file=sys.stderr)
            sys.exit()
        if arg['--pvalue'] == True and 'log10PVALUE' not in arg['--fields']:
            print('Error: is not possible make P-VALUE graphics with your input data', file=sys.stderr)
            sys.exit()
        if arg['--allele-freq-1'] == True and 'SNPidx1' not in arg['--fields']:
            print('Error: is not possible make phased Allele Frequency graphics with your input data. Please use -M or -a option', file=sys.stderr)
            sys.exit()
        if arg['--allele-freq-2'] == True and 'SNPidx2' not in arg['--fields']:
            print('Error: is not possible make Allele Frequency graphics for this pool. Please use -R, -M or -a option', file=sys.stderr)
            sys.exit()
        if arg['--max-allele-freq2'] == True and 'MAX_SNPidx2' not in arg['--fields']:
            print('Error: is not possible make phased MAX Allele Frequency graphics with your input data. Please use -D or -R or -a option', file=sys.stderr)
            sys.exit()
    if arg['--all'] == False and arg['--pvalue'] == False and arg['--allele-freq-1'] == False and arg['--allele-freq-2']:
        print('Warning: any graphic was selected. Please use -a option to make all possible graphics with your data.', file=sys.stderr)
    return arg


def check_qtl_opts(arg, df):
    arg['titles'] = titles_qtl
    arg['lines'] = lines_qtl

    if arg['--all']:
        if 'log10PVALUE' in arg['--fields']:
            arg['--pvalue'] = True
        else:
            arg['--pvalue'] = False
        if arg['--ref-genotype'] == True:
            arg['--allele-freq-H'] = True
            arg['--allele-freq-L'] = True
            arg['--combine'] = True
        if 'DELTA' in arg['--fields']:
            arg['--delta'] = True
            arg['--ci95'] = True
        if len(arg['--chromosomes']) > 1:
            arg['--multi-chrom'] = True
        if 'ED' in arg['--fields']:
            arg['--euclidean-distance'] = True
            for ch in arg['--chromosomes']:
                d = df[df['#CHROM'] == ch]
                if len(d) < RANG:
                    arg['--euclidean-distance'] = False
                    print(f'Warning: there is not enough markers in {ch} to calculate ED100', file=sys.stderr)
        if 'G' in arg['--fields']:
            arg['--g-statistic'] = True
        if arg['--pvalue'] and arg['--delta'] and arg['--euclidean-distance'] and arg['--g-statistic']:
            arg['--qtl-seq'] = True
    else:
        if arg['--pvalue'] == True and 'log10PVALUE' not in arg['--fields']:
            print('Error: is not possible make P-VALUE graphics with your input data', file=sys.stderr)
            sys.exit()
        if arg['--delta'] == True and 'DELTA' not in arg['--fields']:
            print('Error: is not possible make phased SNP-idx graphics with your input data. Please use -a.', file=sys.stderr)
            sys.exit()
        if arg['--ci95'] == True and 'DELTA' not in arg['--fields']:
            print('Error: is not possible to calculate confidense interval without DELTA field. Please use -a.', file=sys.stderr)
            sys.exit()
        if arg['--euclidean-distance'] == True and 'ED' not in arg['--fields']:
            print('Error: is not possible make Euclidean Distance graphics with your input data', file=sys.stderr)
            sys.exit()
        if arg['--g-statistic'] == True and 'G' not in arg['--fields']:
            print('Error: is not possible make G-statistic graphics with your input data', file=sys.stderr)
            sys.exit()
        if (arg['--ref-genotype'] == False) and (arg['--allele-freq-H'] == True or arg['--allele-freq-L'] == True):
            print('Error: is not possible make phased SNP-idx graphics with your input data. Please use -a or -p.', file=sys.stderr)
            sys.exit()
        if arg['--combine'] == True and arg['--ref-genotype'] == False:
            print('Error: is not possible to make combined graphics with your input data. Please use -a.', file=sys.stderr)
            sys.exit()
    if arg['--all'] == False and arg['--pvalue'] == False and arg['--delta'] == False and arg['--allele-freq-H'] == False and arg['--allele-freq-L'] == False and arg['--qtl-seq'] == False and arg['--g-statistic'] == False and arg['--euclidean-distance'] == False:
        print('Warning: any graphic was selected. Please use -a option to make all possible graphics with your data.', file=sys.stderr)

    return arg


def check_save(arg, file_name):
    typ = arg['--output-type']
    if os.path.isfile(arg['--outdir']+file_name):

        expand = 0
        while True:
            expand += 1
            nw_file_name = file_name.split(typ)[0] + '_' + str(expand) + typ
            if os.path.isfile(arg['--outdir']+nw_file_name):
                continue
            else:
                file_name = nw_file_name
                return file_name
    else:
        return file_name


def pvalor2(a, b, c, d):
  # a,b,c,d = row[0], row[1],row[2],row[3]
  table = np.array([[a, b], [c, d]])
  oddsr, pvalue = fisher_exact(table, alternative='two-sided')
  # print(pvalue)
  if pvalue == 0:
    return 0, -np.inf
  else:
    return pvalue, log(pvalue)/log(10)


def expected_pvalue(row):
    a = int(row[0])
    b = int(row[1])
    c = int(row[2])
    d = int(row[3])
    fila1 = a + b
    fila2 = c + d
    col1 = a + c
    col2 = b + d
    tot = col1 + col2
    a_ex = round((fila1*col1)/tot)
    b_ex = round((fila1*col2)/tot)
    c_ex = round((fila2*col1)/tot)
    d_ex = round((fila2*col2)/tot)
    e_pvalue = pvalor(a_ex, b_ex, c_ex, d_ex)
    e_logpvalue = -(log(e_pvalue)/log(10))
    return e_logpvalue


def grouped_by(df, arg):
    w = arg['--window']
    chrom = arg['--chromosomes']
    spacer = arg['spacer']
    fsal = arg['fsal']
    n_head = spacer.join(arg['--fields'])+'\n'
    write_line(n_head, arg['fsal'])
    for ch in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[ch]]
        if arg['--window'] != False:
            d.index = np.arange(len(d))
            d_list = [i for i in range(len(d))]
            g_list = [d_list[i:i+w] for i in range(0, len(d_list), w)]  # indexes for each group
        else:
            D = arg['--distance']
            maxPos = arg['--contigs'][chrom[ch]]
            intRanges = [x for x in range(0, maxPos,D)]
            if intRanges[-1] < maxPos:
                intRanges.append(maxPos)
            intLabels = [f'{start}-{end-1}' for start, end in zip(intRanges, intRanges[1:])]
            d['Interval'] = pd.cut(d['POS'], bins=intRanges, labels=intLabels)
            grouped_indices = d.groupby('Interval').apply(lambda group: group.index.tolist())
            g_list = [indices for indices in grouped_indices if indices]

        for g in g_list:  # for each group in group list
            res = list()
            if 'DPdom_2' not in arg['--fields'] and 'DPrec_2' not in arg['--fields']:
                rAt, rBt, rPost, rTt = 0, 0, 0, 0
            else:
                rAt, rBt, rCt, rDt, rPost, rTt = 0, 0, 0, 0, 0, 0
            for i in g:  # for each index in group
                if 'DPdom_2' not in arg['--fields'] and 'DPrec_2' not in arg['--fields']:
                    rA = d.loc[i].DPdom_1
                    rB = d.loc[i].DPrec_1
                    rT = rA + rB
                    rPos = d.loc[i].POS*rT
                    rTt += rT
                    rAt += rA
                    rBt += rB
                    rPost += rPos
                else:
                    rA = d.loc[i].DPdom_1
                    rB = d.loc[i].DPrec_1
                    rC = d.loc[i].DPdom_2
                    rD = d.loc[i].DPrec_2
                    rT = rA + rB + rC + rD
                    rPos = d.loc[i].POS*rT
                    rTt += rT
                    rAt += rA
                    rBt += rB
                    rCt += rC
                    rDt += rD
                    rPost += rPos
            rPost = round(rPost/(rTt))
            res += [chrom[ch], rPost, rAt, rBt]
            if 'DPdom_2' in arg['--fields'] and 'DPrec_2' in arg['--fields']:
                res += [rCt, rDt]
                if 'SNPidx1' in arg['--fields'] and 'SNPidx2' in arg['--fields']:
                    rSNPidx1 = rBt/(rAt+rBt)
                    rSNPidx2 = rDt/(rCt+rDt)
                    res += [rSNPidx1, rSNPidx2]
                if 'MAX_SNPidx2' in arg['--fields']:
                    rMAX_SNPidx2 = (max(rDt, rCt))/(rDt+rCt)
                    rSNPidx1 = rBt/(rAt+rBt)
                    rSNPidx2 = rDt/(rCt+rDt)
                    fisher = LogFisher(rAt, rBt, rCt, rDt)
                    boost = 1/(sys.float_info.min + abs(1 - 1/max(rSNPidx2, 1-rSNPidx2)))
                    res += [rMAX_SNPidx2, fisher, boost]
                if 'DELTA' in arg['--fields']:
                    rDELTA = rSNPidx2 - rSNPidx1
                    res += [rDELTA]
                if 'G' in arg['--fields'] and 'ED' in arg['--fields']:
                    v1 = (rCt/(rCt + rDt), rDt/(rCt + rDt))
                    v2 = (rAt/(rAt + rBt), rBt/(rAt + rBt))
                    ed = distance.euclidean(v1,v2)
                    g = Gstatic(rAt, rBt, rCt, rDt)
                    res += [ed,g]
                pval, log10pval = pvalor2(rAt, rBt, rCt, rDt)
                res += [pval, log10pval]
            else:
                if 'MAX_SNPidx2' in arg['--fields']:
                    rMAX_SNPidx2 = (max(rAt, rBt))/(rAt+rBt)
                    boost = 1/(sys.float_info.min + abs(1 - 1 /
                               max(rMAX_SNPidx2, 1-rMAX_SNPidx2)))
                    res += [rMAX_SNPidx2, boost]
                elif 'SNPidx1' in arg['--fields']:
                    rSNPidx1 = rBt/(rAt+rBt)
                    boost = 1/(sys.float_info.min +
                               abs(1 - 1/max(rSNPidx1, 1-rSNPidx1)))
                    res += [rSNPidx1, boost]

            n_line = spacer.join(str(field) for field in res) + '\n'
            write_line(n_line, fsal)


def get_ED100_4(df, arg, rang):
    chrom = [ch for ch in arg['contigs'].keys()]
    #chrom = df['#CHROM'].unique()
    ED100 = np.array([])
    for ch in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[ch]]
        if len(d) < RANG:
            continue
        d.index = np.arange(len(d))
        n = len(d)
        ED100c = np.zeros(n)
        for i in range(0, int(rang/2)):
            ED100c[i] = d['ED'].iloc[0:rang].sum()
        for i in range(int(rang/2), n-int(rang/2)):
            ED100c[i] = d['ED'].iloc[i-int(rang/2):i+int(rang/2)].sum()
        for i in range(n-int(rang/2), n):
            ED100c[i] = d['ED'].iloc[n-int(rang/2):n].sum()
        ED100 = np.concatenate([ED100, ED100c])
    df['ED100_4'] = pd.Series(ED100**4)/(10**8)
    return df


def plot_ED(df, arg):
    chrom = arg['--chromosomes']
    typ = arg['--output-type']
    max_y = max(df['ED100_4'])
    t,rt,_=arg['titles'][10]
    for ch in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[ch]]
        d.index = np.arange(len(d))
        max_x = int(arg['contigs'][chrom[ch]])
        x = d[['POS']]
        y = d[['ED']]
        fig, ax = plt.subplots(figsize=(10, 4.2))
        ax.scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax.set(xlim=(0, max_x), ylim=(0, 1.5))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=8)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.set_xlabel(xlabel='Chromosomal position (pb)', fontsize=15)
        ax.set_ylabel(ylabel='Euclidean distance', fontsize=12, rotation=90, labelpad=15)
        ax.set_yticks(ticks=[0, 0.5, 1, 1.5])
        ax.set_yticklabels(labels=[0, 0.5, 1, 1.5], fontsize=8)
        ax.spines['top'].set_visible(False)
        plot_ED100_4(d, arg, ax, max_x, max_y)

        rtch = rt.format(chrom[ch])
        filename = rtch + typ
        filename = check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()

        cap = list()
        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t)
            cap.append(arg['lines'][9].format(arg['color_names']['dots']))
            cap.append(arg['lines'][10].format(arg['color_names']['mvg']))
            write_caption(f,cap,arg)

def plot_G(df, arg):
    chrom = arg['--chromosomes']
    typ = arg['--output-type']
    min_y, max_y =min(df['G']), max(df['G'])
    t,rt,_ = arg['titles'][12]
    for ch in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[ch]]
        max_x = int(arg['contigs'][chrom[ch]])
        x = d[['POS']]
        y = d[['G']]
        fig, ax = plt.subplots(figsize=(10, 4.2))
        ax.scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax.set(xlim=(0, max_x), ylim=(min_y, max_y))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=8)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.set_xlabel(xlabel='Chromosomal position (pb)', fontsize=15)
        ax.set_ylabel(ylabel='G-statistic', fontsize=12, rotation=90, labelpad=15)
        #ax.set_yticks(ticks=[0, 0.5, 1, 1.5])
        #ax.set_yticklabels(labels=[0, 0.5, 1, 1.5], fontsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax, 'G')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[ch], arg, ax, 'G')
        rtch = rt.format(chrom[ch])
        filename = rtch + typ
        filename = check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
        cap = list()

        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t)
            cap.append(arg['lines'][11].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][12].format(arg['color_names']['mvg'], str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][-13].format(arg['color_names']['mvg'], str(arg['--distance-avg'])))
            write_caption(f,cap,arg)

def pval_multi_graph(df, arg):
    t, rt = arg['titles'][1]
    if arg['--captions']:
        f = create_caption(arg, rt)
    typ = arg['--output-type']
    min_y = min(df['log10PVALUE'])*1.05
    labs_list = list()
    f_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig, ax = plt.subplots(1, len(chrom), figsize=(len(chrom)*(3.3333333335), 2))
        if len(chrom) == 1:
            fig, ax = plt.subplots(1,2,figsize=(2*(3), 2))
        for i in range(len(chrom)):
            d = df[df['#CHROM'] == chrom[i]]
            max_x = int(arg['contigs'][chrom[i]])
            x = d[['POS']]
            y = d[['log10PVALUE']]
            ax[i].scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
            ax[i].set(xlim=(0, max_x), ylim=(min_y, 0))
            ax[i].set_xticks([])
            ax[i].set_xlabel(xlabel='chr {}'.format(chrom[i]), fontsize=8)
            labs_list.append('Chromosome {}'.format(chrom[i]))
            ax[i].set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)',fontsize=12,labelpad=15)
            ax[i].tick_params(axis='y', which='major', labelsize=8)
            if arg['--moving-avg'] != False:
                    plot_avg(d, arg, ax[i], 'log10PVALUE')
            if arg['--distance-avg'] != False:
                    plotDistanceAvg(d, chrom[i], arg, ax[i], 'log10PVALUE')
            if arg['--bonferroni']:
                threshold = log(0.05/len(df.axes[0]))/log(10)
                max_x_ch = (max(d['POS'])/max_x)
                ax[i].axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
            if len(chrom) == 1:
                ax[1].remove()
            if i > 0:
                ax[i].axes.get_yaxis().set_visible(False)
        fig.subplots_adjust(wspace=0)
        filename = rt + typ

        filename = check_save(arg, filename)
        f_name.append(filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
    cap = list()
    if arg['--captions']:
        cap.append(', '.join(f_name))
        cap.append(t.format(', '.join(labs_list)))
        cap.append(arg['lines'][0].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][1].format(arg['color_names']['log10PVALUE'], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][-7].format(arg['color_names']['log10PVALUE'], str(arg['--distance-avg'])))
        if arg['--bonferroni']:
            cap.append(arg['lines'][2].format(str(arg['n_markers'])))
        write_caption(f, cap,arg)

def pval_manhattan_plot(df, arg):
    t, rt = arg['titles'][-1]
    if arg['--captions']:
        f = create_caption(arg, rt)
    typ = arg['--output-type']
    max_y = max(-df['log10PVALUE'])*1.05
    labs_list = list()
    f_name = list()
    am = 2.5
    chrom = arg['--chromosomes']
    if len(arg['--chromosomes'])*am <= 18:
        af = len(arg['--chromosomes'])*am
    else:
        af = 18
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, 3))
    for i in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[i]]
        max_x = int(arg['contigs'][chrom[i]])
        x = d[['POS']]
        y = -d[['log10PVALUE']]
        ax[i].scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[i].set(xlim=(0, max_x), ylim=(0, max_y))
        ax[i].set_xticks([])
        ax[i].set_xlabel(xlabel=f'{chrom[i]}', fontsize=12)
        labs_list.append('Chromosome {}'.format(chrom[i]))
        ax[i].set_ylabel(ylabel='-log'+r'$_{10}$'+'(p-value)',fontsize=12,labelpad=15)
        ax[i].tick_params(axis='y', which='major', labelsize=8)
        if arg['--moving-avg'] != False:
            d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
            c_=arg['--palette']['log10PVALUE']
            d['prody']=d['log10PVALUE'] * d['DP']
            d['avgy']=(d['prody'].rolling(arg['--moving-avg']).sum() / \
                d['DP'].rolling(arg['--moving-avg']).sum())
            d['prodx']=d['POS'] * d['DP']
            d['avgx']=(d['prodx'].rolling(arg['--moving-avg']).sum() / \
            d['DP'].rolling(arg['--moving-avg']).sum())
            ax[i].plot(d['avgx'], -d['avgy'], c=c_, lw=2)#lw=0.75
        if arg['--distance-avg'] != False:
            d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
            c_=arg['--palette']['log10PVALUE']
            d['prody']=d['log10PVALUE'] * d['DP']
            d['prodx'] = d['DP']*d['POS']
            interval_ranges = [x for x in range(0, max_x, arg['--distance-avg'])]
            if interval_ranges[-1] < max_x:
                interval_ranges.append(max_x)
            interval_labels = [f'{start}-{end-1}' for start, end in zip(interval_ranges, interval_ranges[1:])]
            d['Interval'] = pd.cut(d['POS'], bins=interval_ranges, labels=interval_labels)
            res = d.groupby('Interval').agg({'prodx': sum, 'prody': sum, 'DP': sum})
            res['POSx'] = (res['prodx']/res['DP'])
            res['VALy'] = (res['prody']/res['DP'])
            res = res.dropna()
            ax[i].plot(res['POSx'], -res['VALy'], c=c_, lw=2)
        if arg['--bonferroni']:
            threshold = -log(0.05/len(df.axes[0]))/log(10)
            max_x_ch = (max(d['POS'])/max_x)
            ax[i].axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        if len(chrom) == 1:
            ax[1].remove()
        if i > 0:
            ax[i].axes.get_yaxis().set_visible(False)
    fig.subplots_adjust(wspace=0)
    filename = rt + typ
    filename = check_save(arg, filename)
    f_name.append(filename)
    plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
    plt.close()
    cap = list()
    if arg['--captions']:
        cap.append(', '.join(f_name))
        cap.append(t.format(', '.join(labs_list)))
        cap.append(arg['lines'][-5].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][-6].format(arg['color_names']['log10PVALUE'], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][-8].format(arg['color_names']['log10PVALUE'], str(arg['--distance-avg'])))
        if arg['--bonferroni']:
            cap.append(arg['lines'][2].format(str(arg['n_markers'])))
        write_caption(f, cap,arg) 

def pval_mono_graph(df, arg):
    chrom=arg['--chromosomes']
    typ=arg['--output-type']
    
    min_y=min(df['log10PVALUE'])*1.05
    
    for i in range(len(chrom)):
        d=df[df['#CHROM'] == chrom[i]]
        max_x= int(arg['contigs'][chrom[i]])
        #DMAX = d[d.DELTA == d.DELTA.max()]
        #DMAX = DMAX[DMAX.log10PVALUE == DMAX.log10PVALUE.min()]
        #print('DELTA__MAX', DMAX)
        #DMIN = d[d.DELTA == d.DELTA.min()]
        #print('DELTA__MIN', DMIN)
        t, rt = arg['titles'][2]
        cap = list()
        # TODO - Save the moving average values
        # mov_avg = pd.DataFrame(d, columns=['mediamovil','mediamovilx'])
        # mov_avg.dropna()
        # mov_avg.to_csv(path_or_buf='./mov_avg_pvalue.txt',sep='\t', header=True)
        x=d[['POS']]
        y=d[['log10PVALUE']]
        fig, ax=plt.subplots(figsize=(10, 4.2))
        ax.scatter(x, y, s=0.5, color=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax.set(xlim=(0, max_x), ylim=(min_y, 0))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=8)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)
        ax.set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)', fontsize=12, labelpad=15)
        ax.tick_params(axis='y', which='major', labelsize=8)
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax, 'log10PVALUE')
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax, 'log10PVALUE')
        if arg['--bonferroni']:
            threshold=log(0.05/len(df.axes[0]))/log(10)
            max_x_ch=(max(d['POS'])/max_x)
            ax.axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        rtch = rt.format(chrom[i])

        filename=rtch + typ
        filename=check_save(arg, filename)

        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()

        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t.format(chrom[i]))
            cap.append(arg['lines'][0].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][1].format(arg['color_names']['log10PVALUE'], str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][-7].format(arg['color_names']['log10PVALUE'], str(arg['--distance-avg'])))
            if arg['--bonferroni']:
                cap.append(arg['lines'][2].format(str(arg['n_markers'])))
            write_caption(f, cap,arg)



def AF1_AF2_mono_graph(df, arg):
    chrom=arg['--chromosomes']
    typ=arg['--output-type']

    for i in range(len(chrom)):
        d=df[df['#CHROM'] == chrom[i]]
        t1,rt1,_=arg['titles'][3]
        t2,rt2,_=arg['titles'][4]
        max_x= int(arg['contigs'][chrom[i]])
        x=d[['POS']]
        y1=d[['SNPidx1']]
        y2=d[['SNPidx2']]
        # AF1
        fig, ax=plt.subplots(figsize=(10, 4.2))
        ax.scatter(x, y1, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax.set(xlim=(0, max_x), ylim=(0, 1))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=8)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.set_xlabel(xlabel='Chromosomal position (bp)'.format(chrom[i]), fontsize=15)
        ax.set_ylabel(ylabel='Allele Frequency', fontsize=12, rotation=90, loc='center', labelpad=15)
        ax.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
        ax.set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=8)
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax, 'SNPidx2')
            plot_avg(d, arg, ax, 'SNPidx1')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax, 'SNPidx2')
            plotDistanceAvg(d, chrom[i], arg, ax, 'SNPidx1')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        rtch1 = rt1.format(chrom[i])
        filename=rtch1 + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
        cap = list()
        if arg['--captions']:
            f = create_caption(arg, rtch1)
            cap.append(filename)
            cap.append(t1.format(chrom[i]))
            cap.append(arg['lines'][-2].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][3].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])))
                cap.append(arg['lines'][4].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][-9].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])))
                cap.append(arg['lines'][-10].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
            write_caption(f,cap,arg) 
        # AF2
        fig, ax=plt.subplots(figsize=(10, 4.2))
        ax.scatter(x, y2, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax.set(xlim=(0, max_x), ylim=(0, 1))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=8)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)
        ax.set_ylabel(ylabel='Allele frequency', fontsize=12, rotation=90, loc='center')
        ax.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
        ax.set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=8)
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax, 'SNPidx1')
            plot_avg(d, arg, ax, 'SNPidx2')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax, 'SNPidx1')
            plotDistanceAvg(d, chrom[i], arg, ax, 'SNPidx2')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        rtch2 = rt2.format(chrom[i])
        filename=rtch2 + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
        cap = list()
        if arg['--captions']:
            f = create_caption(arg, rtch2)
            cap.append(filename)
            cap.append(t2.format(chrom[i]))
            cap.append(arg['lines'][-1].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][4].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
                cap.append(arg['lines'][3].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][-9].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
                cap.append(arg['lines'][-10].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])))
            write_caption(f,cap,arg)


def AF_mono_graph(df, arg, g_type):
    chrom=arg['--chromosomes']
    typ=arg['--output-type']
    ticks_y=[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    ylab = 'Allele Frequency'
    boost = False
    for i in range(len(chrom)):
        if g_type == 'SNPidx1':
            key2 = -2 
            if arg['--moving-avg'] != False:
                key = 3
            if arg['--distance-avg'] != False:
                key = -9
            t,_,rt=arg['titles'][3]
        if g_type == 'SNPidx2':
            key2 = -1
            if arg['--moving-avg'] != False:
                key = 4
            if arg['--distance-avg'] != False:
                key = -10
            t,_,rt=arg['titles'][4]
        if g_type == 'MAX_SNPidx2':
            ylab = 'Maximum Allele Frequency'
            key2 = -3
            if arg['--moving-avg'] != False:
                key = 5
            if arg['--distance-avg'] != False:
                key = -11
            t,rt=arg['titles'][5]
            ticks_y=[0.5, 0.75, 1]
            lim_y=(0.5, 1)
        if g_type == 'DELTA':
            if arg['--ref-genotype'] == True:
                ylab = '$\Delta$'+' (SNP-index)'
                key2 = -3
                if arg['--moving-avg'] != False:
                    key = 5
                if arg['--distance-avg'] != False:
                    key = -11
                t,rt = arg['titles'][5]
                ticks_y=[-1, 0, 1]
                lim_y=(-1.2, 1.2)
            else:
                ylab = '|$\Delta$'+' (SNP-index)|'
                key2 = -4
                if arg['--moving-avg'] != False:
                    key = 5
                if arg['--distance-avg'] != False:
                    key = -11
                t,rt = arg['titles'][5]
                ticks_y=[0, 0.25, 0.5, 0.75, 1]
                lim_y=(0, 1.2)
        d=df[df['#CHROM'] == chrom[i]]
        max_x=int(arg['contigs'][chrom[i]])
        x=d[['POS']]
        y=d[[g_type]]
        fig, ax=plt.subplots(figsize=(10, 4.2))
        ax.scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax.set(xlim=(0, max_x), ylim=lim_y)
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=8)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)
        ax.set_ylabel(ylabel=ylab, fontsize=12, rotation=90, labelpad=15)
        ax.set_yticks(ticks=ticks_y)
        ax.set_yticklabels(labels=ticks_y, fontsize=8)
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax, g_type)
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax, g_type)
        if 'mbsplot' in arg.keys():
            if arg['--boost'] != False:
                keyB = 6
                if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                    plot_boost(d, arg, ax, max_x)
                    boost = True
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plot_boost(d, arg, ax, max_x)
                    boost = True
            if arg['--distance-boost'] != False:
                keyB = -12
                if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                    plotDistanceBoost(d, arg, ax, max_x, chrom[i])
                    boost = True
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plotDistanceBoost(d, arg, ax, max_x, chrom[i])
                    boost = True
        if 'qtlplot' in arg.keys():
            if arg['--ci95'] and g_type == 'DELTA':
                if arg['--moving-avg'] != False:
                    calc_ci(d, arg, ax)
                if arg['--distance-avg'] != False:
                    distanceCI(d, arg, ax, chrom[i])
                ax.axhline(y=0, color='black', linestyle='dashed', linewidth=0.75)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        rtch = rt.format(chrom[i])
        filename=rtch + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
        cap = list()
        if arg['--captions']:
            f = create_caption(arg,rtch)
            cap.append(filename)
            cap.append(t.format(chrom[i]))
            cap.append(arg['lines'][key2].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][key].format(arg['color_names'][g_type], str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][key].format(arg['color_names'][g_type], str(arg['--distance-avg'])))
            if 'mbsplot' in arg.keys():
                if boost == True:
                    if arg['--boost'] != False:
                        cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--boost'])))
                    if arg['--distance-boost'] != False:
                        cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--distance-boost'])))
            if 'qtlplot' in arg.keys():
                if arg['--ci95'] and g_type == 'DELTA':
                    cap.append(arg['lines'][7].format(arg['color_names']['ci'], str(arg['n_markers'])))
            write_caption(f, cap,arg)

def pval_multi_Vertical_graph(df, arg):
    typ=arg['--output-type']
    min_y = min(df['log10PVALUE'])*1.05
    max_x = arg['max_lenght']
    t,rt = arg['titles'][6]
    labs_list = list()
    f_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig,ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['log10PVALUE']]
            ax[i].scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
            ax[i].set(xlim=(0, max_x), ylim=(min_y, 0))
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)', fontsize=12, rotation=90, labelpad=15)
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=16, rotation=0, x=-0.12, y=0.85)
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].tick_params(axis='y', which='major', labelsize=8)
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], 'log10PVALUE')
            if arg['--distance-avg'] != False:
                plotDistanceAvg(d, chrom[i], arg, ax[i], 'log10PVALUE')
            if arg['--bonferroni']:
                threshold=log(0.05/len(df.axes[0]))/log(10)
                max_x_ch=(max(d['POS'])/max_x)
                ax[i].axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            if len(chrom) == 1:
                ax[1].remove()
            if i == len(chrom)-1 or len(chrom) == 1:
                ax[i].tick_params(labelbottom=True)
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=8, color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)


        fig.subplots_adjust(hspace=0.1)
        filename= rt + typ
        filename=check_save(arg, filename)
        f_name.append(filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        cap = cap + labs_list
        cap.append(arg['lines'][0].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][1].format(arg['color_names']['log10PVALUE'],str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
                cap.append(arg['lines'][-7].format(arg['color_names']['log10PVALUE'], str(arg['--distance-avg'])))
        if arg['--bonferroni']:
            cap.append(arg['lines'][2].format(str(arg['n_markers'])))
        write_caption(f,cap,arg)


def ED_multi_Vertical_graph(df, arg):
    typ=arg['--output-type']
    max_y = max(df['ED100_4'])
    max_x = arg['max_lenght']
    t,_,rt = arg['titles'][10]
    labs_list = list()
    f_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig,ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['ED']]
            ax[i].scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
            ax[i].set(xlim=(0, max_x), ylim=(0, 1.5))
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='Euclidean distance', fontsize=12, rotation=90, labelpad=15)
            ax[i].set_yticks(ticks=[0,0.5,1,1.5])
            ax[i].set_yticklabels(labels=[0,0.5,1,1.5], fontsize=8)
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=16, rotation=0, x=-0.12, y=0.85)
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].tick_params(axis='y', which='major', labelsize=8)
            ax[i].spines['top'].set_visible(False)
            plot_ED100_4(d, arg, ax[i], max_x, max_y)
            if len(chrom) == 1:
                ax[1].remove()
            if i == len(chrom)-1 or len(chrom) == 1:
                ax[i].tick_params(labelbottom=True)
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=8, color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)

        fig.subplots_adjust(hspace=0.2)
        filename= rt + typ
        filename=check_save(arg, filename)
        f_name.append(filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        cap = cap + labs_list
        cap.append(arg['lines'][9].format(arg['color_names']['dots']))
        cap.append(arg['lines'][10].format(arg['color_names']['mvg']))
        write_caption(f,cap,arg)

def G_multi_Vertical_graph(df, arg):
    typ=arg['--output-type']
    min_y, max_y = min(df['G']), max(df['G'])
    max_x =arg['max_lenght']
    t,_,rt = arg['titles'][12]
    labs_list = list()
    f_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig,ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['G']]
            ax[i].scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
            ax[i].set(xlim=(0, max_x), ylim=(min_y, max_y))
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='G-statistic', fontsize=12, rotation=90, labelpad=15)
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=16, rotation=0, x=-0.12, y=0.85)
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].tick_params(axis='y', which='major', labelsize=8)
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], 'G')
            if arg['--distance-avg'] != False:
                plotDistanceAvg(d, chrom[i], arg, ax[i], 'G')
            if len(chrom) == 1:
                ax[1].remove()
            if i == len(chrom)-1 or len(chrom) == 1:
                ax[i].tick_params(labelbottom=True)
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=8, color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)

        fig.subplots_adjust(hspace=0.1)
        filename= rt + typ
        filename=check_save(arg, filename)
        f_name.append(filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        cap = cap + labs_list
        cap.append(arg['lines'][11].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][12].format(arg['color_names']['mvg'], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
                cap.append(arg['lines'][-13].format(arg['color_names']['mvg'], str(arg['--distance-avg'])))
        write_caption(f,cap,arg)

def AF_multi_Vertical_graph(df, arg, g_type):
    typ=arg['--output-type']
    max_x=arg['max_lenght']
    ticks_y=[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    ylab = 'Allele Frequency'
    boost = False
    labs_list = list()
    f_name = list()
    if g_type == 'SNPidx1':
        key2 = -2
        if arg['--moving-avg'] != False:
            key = 3
        if arg['--distance-avg'] != False:
            key = -9
        t,_,rt=arg['titles'][7]
    if g_type == 'SNPidx2':
        key2 = -1
        if arg['--moving-avg'] != False:
            key = 4
        if arg['--distance-avg'] != False:
            key = -10
        t,_,rt=arg['titles'][8]
    if g_type == 'MAX_SNPidx2':
        ylab = 'Maximum Allele Frequency'
        key2 = -3
        if arg['--moving-avg'] != False:
            key = 5
        if arg['--distance-avg'] != False:
            key = -11
        t,rt=arg['titles'][9]
        ticks_y=[0.5, 0.75, 1]
        lim_y=(0.5, 1)
    if g_type == 'DELTA':
        if arg['--ref-genotype'] == True:
            ylab = '$\Delta$' + '(SNP-index)'
            key2 = -3
            if arg['--moving-avg'] != False:
                key = 5
            if arg['--distance-avg'] != False:
                key = -11
            t,rt = arg['titles'][9]
            ticks_y = [-1,0,1]
            lim_y = (-1.2, 1.2)
        else:
            ylab = '|$\Delta$' + '(SNP-index)|'
            key2 = -4
            if arg['--moving-avg'] != False:
                key = 5
            if arg['--distance-avg'] != False:
                key = -11
            t,rt = arg['titles'][9]
            ticks_y = [0,0.25,0.5,0.75,1]
            lim_y = (0, 1.2)

    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig, ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10,2.5*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[[g_type]]

            ax[i].scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
            ax[i].set(xlim=(0, max_x), ylim=lim_y)
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel=ylab, fontsize=12, rotation=90, labelpad=15)
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=16, rotation=0, x = -0.12, y=0.85)
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].set_yticks(ticks=ticks_y)
            ax[i].set_yticklabels(labels=ticks_y, fontsize=8)
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], g_type)
            if arg['--distance-avg'] != False:
                plotDistanceAvg(d, chrom[i], arg, ax[i], g_type)
            if 'mbsplot' in arg.keys():
                if arg['--boost'] != False:
                    keyB = 6
                    if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                        plot_boost(d, arg, ax[i], max_x)
                        boost = True
                    if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                        plot_boost(d, arg, ax[i], max_x)
                        boost = True
                if arg['--distance-boost'] != False:
                    keyB = -12
                    if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                        plotDistanceBoost(d, arg, ax[i], max_x, chrom[i])
                        boost = True
                    if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                        plotDistanceBoost(d, arg, ax[i], max_x, chrom[i])
                        boost = True
            if 'qtlplot' in arg.keys():
                if arg['--ci95'] and g_type == 'DELTA':
                    if arg['--moving-avg'] != False:
                        calc_ci(d, arg, ax[i])
                    if arg['--distance-avg'] != False:
                        distanceCI(d, arg, ax[i], chrom[i])
                    ax[i].axhline(y=0, color='black', linestyle='dashed', linewidth=0.75)

            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            if len(chrom) == 1:
                ax[1].remove()
            if i == len(chrom)-1 or len(chrom) == 1:
                ax[i].tick_params(labelbottom=True)
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=8)
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)

        fig.subplots_adjust(hspace=0.1)
        filename= rt + typ
        filename=check_save(arg, filename)
        f_name.append(filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        cap = cap + labs_list
        cap.append(arg['lines'][key2].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][key].format(arg['color_names'][g_type], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][key].format(arg['color_names'][g_type], str(arg['--distance-avg'])))
        if 'mbsplot' in arg.keys():
            if boost == True:
                if arg['--boost'] != False:
                    cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--boost'])))
                if arg['--distance-boost'] != False:
                    cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--distance-boost'])))
        if 'qtlplot' in arg.keys() and g_type == 'DELTA':
            if arg['--ci95']:
                cap.append(arg['lines'][7].format(arg['color_names']['ci'], str(arg['n_markers'])))
        write_caption(f, cap,arg)

def AF12_multi_Vertical_graph(df, arg):
    typ=arg['--output-type']
    max_x=arg['max_lenght']
    
    ticks_y=[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    t1,rt1,_=arg['titles'][7]
    t2,rt2,_=arg['titles'][8]
    labs_list = list()
    f1_name = list()
    f2_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig, ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['SNPidx1']]
            #AF1&2
            ax[i].scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
            ax[i].set(xlim=(0, max_x), ylim=lim_y)
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='Allele Frequency', fontsize=12, rotation=90, labelpad=15)
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=16, rotation=0, x = -0.12, y=0.85) 
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))  
            ax[i].set_yticks(ticks=ticks_y)
            ax[i].set_yticklabels(labels=ticks_y, fontsize=8)
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], 'SNPidx2')
                plot_avg(d, arg, ax[i], 'SNPidx1')
            if arg['--distance-avg'] != False:
                plotDistanceAvg(d, chrom[i], arg, ax[i], 'SNPidx2')
                plotDistanceAvg(d, chrom[i], arg, ax[i], 'SNPidx1')
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            if len(chrom) == 1:
                ax[1].remove()
            if i == len(chrom)-1 or len(chrom) == 1:
                ax[i].tick_params(labelbottom=True)
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=8, color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)

        fig.subplots_adjust(hspace=0.1)
        filename=rt1+typ
        filename=check_save(arg, filename)
        f1_name.append(filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt1)
        cap.append(', '.join(f1_name))
        cap.append(t1)
        cap = cap + labs_list
        cap.append(arg['lines'][-2].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][3].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])))
            cap.append(arg['lines'][4].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][-9].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])))
            cap.append(arg['lines'][-10].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
        write_caption(f,cap,arg)
    
    #AF1&2
    f2_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig,ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['SNPidx2']]
            ax[i].scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
            ax[i].set(xlim=(0, max_x), ylim=lim_y)
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='Allele Frequency', fontsize=12, rotation=90, labelpad=15)
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=16, rotation=0, x = -0.12, y=0.85)
            ax[i].set_yticks(ticks=ticks_y)
            ax[i].set_yticklabels(labels=ticks_y, fontsize=8)
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], 'SNPidx1')
                plot_avg(d, arg, ax[i], 'SNPidx2')
            if arg['--distance-avg'] != False:
                plotDistanceAvg(d, chrom[i], arg, ax[i], 'SNPidx1')
                plotDistanceAvg(d, chrom[i], arg, ax[i], 'SNPidx2')
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            if len(chrom) == 1:
                ax[1].remove()
            if i == len(chrom)-1 or len(chrom) == 1:
                ax[i].tick_params(labelbottom=True)
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=8, color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)

        fig.subplots_adjust(hspace=0.1)
        filename=rt2+typ
        filename=check_save(arg, filename)
        f2_name.append(filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()

    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt2)
        cap.append(', '.join(f2_name))
        cap.append(t2)
        cap = cap + labs_list
        cap.append(arg['lines'][-1].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][4].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
            cap.append(arg['lines'][3].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][-9].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
            cap.append(arg['lines'][-10].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])))
        write_caption(f,cap,arg)

def AF1_AF2_pval_mono(df, arg):
    global fields
    chrom=arg['--chromosomes']
    typ=arg['--output-type']
    min_y=min(df['log10PVALUE'])*1.05

    for i in range(len(chrom)):
        fig, ax=plt.subplots(3, 1, figsize=(8.5, 7))#(7,9)paper
        d=df[df['#CHROM'] == chrom[i]]
        max_x = arg['contigs'][chrom[i]]
        t,rt = arg['titles'][10]
        x=d[['POS']]
        y1=d['SNPidx1']
        y2=d['SNPidx2']
        y3=d['log10PVALUE']
        # SNPidx1
        ax[0].scatter(x, y1, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[0].set(xlim=(0, max_x), ylim=(0, 1))
        ax[0].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[0].tick_params(labelbottom=False)
        ax[0].set_ylabel(ylabel='Allele Frequency', fontsize=12, rotation=90, labelpad=15)
        ax[0].set_title('(a)', fontsize=16, rotation=0, x = -0.14, y=0.85)
        ax[0].set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
        ax[0].set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=8)
        if arg['--moving-avg'] != False:
            if arg['--combine']:
                plot_avg(d, arg, ax[0], 'SNPidx2')
                plot_avg(d, arg, ax[0], 'SNPidx1')
            else:
               plot_avg(d, arg, ax[0], 'SNPidx1')
        
        if arg['--distance-avg'] != False:
            if arg['--combine']:
                plotDistanceAvg(d, chrom[i], arg, ax[0], 'SNPidx2')
                plotDistanceAvg(d, chrom[i], arg, ax[0], 'SNPidx1')
            else:
               plotDistanceAvg(d, chrom[i], arg, ax[0], 'SNPidx1')
        ax[0].spines['top'].set_visible(False)
        ax[0].spines['right'].set_visible(False)

        # SNPidx2
        ax[1].scatter(x, y2, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[1].set(xlim=(0, max_x), ylim=(0, 1))
        ax[1].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[1].tick_params(labelbottom=False)
        ax[1].set_ylabel(ylabel='Allele Frequency', fontsize=12, rotation=90, labelpad=15)
        ax[1].set_title('(b)', fontsize=16, rotation=0, x = -0.14, y=0.85)
        ax[1].set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
        ax[1].set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=8)
        if arg['--moving-avg'] != False:
            if arg['--combine']:
                plot_avg(d, arg, ax[1], 'SNPidx1')
                plot_avg(d, arg, ax[1], 'SNPidx2')
            else:
               plot_avg(d, arg, ax[1], 'SNPidx2')
        if arg['--distance-avg'] != False:
            if arg['--combine']:
                plotDistanceAvg(d, chrom[i], arg, ax[1], 'SNPidx1')
                plotDistanceAvg(d, chrom[i], arg, ax[1], 'SNPidx2')
            else:
               plotDistanceAvg(d, chrom[i], arg, ax[1], 'SNPidx2')
        ax[1].spines['top'].set_visible(False)
        ax[1].spines['right'].set_visible(False)

        # log10PVALUE
        ax[2].scatter(x, y3, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[2].set(xlim=(0, max_x), ylim=(min_y, 0))
        ax[2].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[2].set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=8)
        ax[2].set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)',fontsize=12, labelpad=15)
        ax[2].set_title('(c)', fontsize=16, rotation=0, x = -0.14, y=0.85)
        ax[2].tick_params(axis='y', which='major', labelsize=8)
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[2], 'log10PVALUE')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[2], 'log10PVALUE')
        if arg['--bonferroni']:
            threshold=log(0.05/len(df.axes[0]))/log(10)
            max_x_ch=(max(d['POS'])/max_x)
            ax[2].axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        ax[2].xaxis.set_major_formatter(ScalarFormatter())
        ax[2].ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax[2].set_xlabel(xlabel='Chromosomal position (bp)',fontsize=15)

        ax[2].spines['top'].set_visible(False)
        ax[2].spines['right'].set_visible(False)

        # Save file
        fig.subplots_adjust(hspace=0.1)
        rtch = rt.format(chrom[i])
        filename=rtch + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
        cap = list()

        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t.format(chrom[i]))
            cap.append(arg['lines'][7].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][3].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])))
                cap.append(arg['lines'][4].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
                cap.append(arg['lines'][1].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][-9].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])))
                cap.append(arg['lines'][-10].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
                cap.append(arg['lines'][-7].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
            if arg['--bonferroni']:
                cap.append(arg['lines'][2].format(str(arg['n_markers'])))
            write_caption(f,cap,arg)

def qtl_mixed_plot(df, arg):
    global fields
    chrom=arg['--chromosomes']
    typ=arg['--output-type']
    t,rt = arg['titles'][13]
    for i in range(len(chrom)):
        fig, ax=plt.subplots(4, 1, figsize=(8.5, 8.5))#(7,9)paper
        d=df[df['#CHROM'] == chrom[i]]
        max_x = arg['contigs'][chrom[i]]
        x=d[['POS']]
        y1=d['G']
        y2=d['ED']
        y3=d['DELTA']
        y4=d['log10PVALUE']
        # G-statistic
        ax[0].scatter(x, y1, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[0].set(xlim=(0, max_x), ylim=(min(df['G']), max(df['G'])))
        ax[0].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[0].tick_params(labelbottom=False)
        ax[0].set_ylabel(ylabel='G-statistic', fontsize=12, rotation=90)
        ax[0].yaxis.set_label_coords(-0.075,0.5)
        ax[0].tick_params(axis='y', which='major', labelsize=8, size=8)
        ax[0].set_title('(a)', fontsize=16, rotation=0, x = -0.14, y=0.85)
        
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[0], 'G' )
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[0], 'G')
        ax[0].spines['top'].set_visible(False)
        ax[0].spines['right'].set_visible(False)

        # ED
        ax[1].scatter(x, y2, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[1].set(xlim=(0, max_x), ylim=(0, 1.5))
        ax[1].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[1].tick_params(labelbottom=False)
        ax[1].set_ylabel(ylabel='Euclidean distance', fontsize=12, rotation=90)
        ax[1].yaxis.set_label_coords(-0.075,0.5)
        ax[1].set_title('(b)', fontsize=16, rotation=0, x = -0.14, y=0.85)
        ax[1].set_yticks(ticks=[0, 0.5, 1, 1.5])
        ax[1].set_yticklabels(labels=[0, 0.5, 1, 1.5], fontsize=8)
        plot_ED100_4(d, arg, ax[1], max_x, max(df['ED100_4']))
        ax[1].spines['top'].set_visible(False)
        ax[1].spines['right'].set_visible(False)

        # DELTA
        if arg['--ref-genotype'] == True:
            ylab = '$\Delta$'+' (SNP-index)'
            ticks_y=[-1, 0, 1]
            lim_y=(-1.2, 1.2)
        else:
            ylab = '|$\Delta$'+' (SNP-index)|'
            ticks_y=[0, 0.25, 0.5, 0.75, 1]
            lim_y=(0, 1.2)

        ax[2].scatter(x, y3, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[2].set(xlim=(0, max_x), ylim=lim_y)
        ax[2].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[2].tick_params(labelbottom=False)
        ax[2].set_ylabel(ylabel= ylab ,fontsize=12)
        ax[2].yaxis.set_label_coords(-0.075,0.5)
        ax[2].set_title('(c)', fontsize=16, rotation=0, x = -0.14, y=0.85)
        ax[2].set_yticks(ticks=ticks_y)
        ax[2].set_yticklabels(labels=ticks_y, fontsize=8)
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[2], 'DELTA')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[2], 'DELTA')
        if arg['--ci95'] and arg['--moving-avg'] != False:
            if arg['--moving-avg'] != False:
                calc_ci(d, arg, ax[2])
            if arg['--distance-avg'] != False:
                distanceCI(d, arg, ax[2], chrom[i])
            ax[2].axhline(y=0, color = 'black', linestyle='dashed', linewidth=0.75)
        ax[2].spines['top'].set_visible(False)
        ax[2].spines['right'].set_visible(False)

        # log10PVALUE
        ax[3].scatter(x, y4, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[3].set(xlim=(0, max_x), ylim=(min(df['log10PVALUE'])*1.05, 0))
        ax[3].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[3].set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=8)
        ax[3].set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)',fontsize=12)
        ax[3].yaxis.set_label_coords(-0.075,0.5)
        ax[3].set_title('(d)', fontsize=16, rotation=0, x = -0.14, y=0.85)
        ax[3].tick_params(axis='y', which='major', labelsize=8)
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[3], 'log10PVALUE')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[3], 'log10PVALUE')
        if arg['--bonferroni']:
            threshold=log(0.05/len(df.axes[0]))/log(10)
            max_x_ch=(max(d['POS'])/max_x)
            ax[3].axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        ax[3].xaxis.set_major_formatter(ScalarFormatter())
        ax[3].ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax[3].set_xlabel(xlabel='Chromosomal position (bp)',fontsize=14)
        ax[3].spines['top'].set_visible(False)
        ax[3].spines['right'].set_visible(False)

        # Save file
        fig.subplots_adjust(hspace=0.1)
        rtch = rt.format(chrom[i])
        filename=rtch + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
        cap = list()

        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t.format(chrom[i]))
            if arg['--ref-genotype'] == True:
                cap.append(arg['lines'][13])
            else:
               cap.append(arg['lines'][14])
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][8].format(arg['color_names']['mvg'], str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][-12].format(arg['color_names']['mvg'], str(arg['--distance-avg'])))
            if arg['--bonferroni']:
                cap.append(arg['lines'][2].format(str(arg['n_markers'])))
            write_caption(f,cap, arg)

def snp_index_graph(df, arg):
    global fields
    chrom=arg['--chromosomes']
    typ=arg['--output-type']
    t,rt=arg['titles'][11]

    for i in range(len(chrom)):
        fig, ax=plt.subplots(3, 1, figsize=(7.5, 5))
        d=df[df['#CHROM'] == chrom[i]]
        max_x=arg['contigs'][chrom[i]]
        x=d[['POS']]
        y1=d['SNPidx1']
        y2=d['SNPidx2']
        y3=d['DELTA']  # delta

        #SNPidx1
        ax[0].scatter(x, y1, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[0].set(xlim=(0, max_x), ylim=(0, 1))
        ax[0].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[0].tick_params(labelbottom=False)
        ax[0].set_yticks(ticks=[0, 0.5, 1])
        ax[0].set_yticklabels(labels=[0, 0.5, 1], rotation=90, fontsize=8)
        ax[0].set_ylabel('SNP-index', fontsize=12)
        ax[0].yaxis.set_label_coords(-0.06,0.5)
        ax[0].set_title('(a)', fontsize=16, rotation=0, x=-0.13, y=0.85)
        ax[0].axhline(y=0.5, color='black', linestyle='dashed', linewidth=0.75)
        ax[0].spines['top'].set_visible(False)
        ax[0].spines['right'].set_visible(False)

        # SNPidx2
        ax[1].scatter(x, y2, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[1].set(xlim=(0, max_x), ylim=(0, 1))
        ax[1].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[1].tick_params(labelbottom=False)
        ax[1].set_yticks(ticks=[0, 0.5, 1])
        ax[1].set_yticklabels(labels=[0, 0.5, 1], fontsize=8, rotation=90)
        ax[1].set_ylabel('SNP-index', fontsize=12)
        ax[1].yaxis.set_label_coords(-0.06,0.5)
        ax[1].set_title('(b)', fontsize=16, rotation=0, x=-0.13, y=0.85)
        ax[1].axhline(y=0.5, color='black', linestyle='dashed', linewidth=0.75)
        ax[1].spines['top'].set_visible(False)
        ax[1].spines['right'].set_visible(False)

        # D-SNPidx
        if arg['--ci95']:
            lim_y = (-1.2, 1.2)
        else:
            lim_y = (-1, 1)
        ax[2].scatter(x, y3, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[2].set(xlim=(0, max_x), ylim=lim_y)
        ax[2].set_yticks(ticks=[-1, 0, 1])
        ax[2].tick_params(labelbottom=True)
        ax[2].set_xticks(ticks=np.arange(0,max_x,5e6))
        ax[2].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=8)
        ax[2].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)
        ax[2].xaxis.set_major_formatter(ScalarFormatter())
        ax[2].ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax[2].set_yticklabels(labels=[-1, 0, 1], rotation=90,fontsize=8)
        ax[2].set_ylabel('$\Delta$'+' (SNP-index)',fontsize=12)
        ax[2].yaxis.set_label_coords(-0.06,0.5)
        ax[2].set_title('(c)', fontsize=16, rotation=0, x=-0.13, y=0.85)
        ax[2].axhline(y=0, color='black', linestyle='dashed', linewidth=0.75)
        ax[2].spines['top'].set_visible(False)
        ax[2].spines['right'].set_visible(False)

        if arg['--moving-avg'] != False:
            plot_avg_qtl_SNPidx(d, arg, ax)
        if arg['--distance-avg'] != False:
            distanceSNPidx(d, arg, ax, chrom[i])
        if arg['--ci95']:
            if arg['--moving-avg'] != False:
                calc_ci(d, arg, ax[2])
            if arg['--distance-avg'] != False:
                distanceCI(d, arg, ax[2], chrom[i])
        # Save file
        fig.subplots_adjust(hspace=0.1)
        rtch = rt.format(chrom[i])
        filename=rtch + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()
        cap = list()
        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t.format(chrom[i]))
            cap.append(arg['lines'][6])
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][8].format(arg['color_names']['mvg'], str(arg['--moving-avg'])))
                if arg['--ci95']:
                    cap.append(arg['lines'][7].format(arg['color_names']['ci'], str(arg['n_markers'])))
            write_caption(f,cap, arg)

def calc_ci(d, arg, ax):
    z95=abs(st.norm.ppf(.025/arg['n_markers']))
    d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
    d['DP1']=d['DPdom_1']+d['DPrec_1']
    d['DP2']=d['DPdom_2']+d['DPrec_2']
    d['productx']=d['POS'] * d['DP']
    d['mediamovilx']=(d['productx'].rolling(arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())

    d['CI95u']=d['DELTA'] + z95 * np.sqrt(d['SNPidx2']*(1-d['SNPidx2'])/d['DP2'] + d['SNPidx1']*(1-d['SNPidx1'])/d['DP1'])
    d['CI95l']=d['DELTA'] - z95 *np.sqrt(d['SNPidx2']*(1-d['SNPidx2'])/d['DP2'] + d['SNPidx1']*(1-d['SNPidx1'])/d['DP1'])

    d['pCIU95']=d['CI95u'] * d['DP']
    d['pCIL95']=d['CI95l'] * d['DP']

    d['Avg95u']=(d['pCIU95'].rolling(arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())
    d['Avg95l']=(d['pCIL95'].rolling(arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())
    
    c_ = arg['--palette']['ci']
    ax.plot(d['mediamovilx'], d['Avg95u'],color=c_, alpha=0.75, lw=0.4)
    ax.plot(d['mediamovilx'], d['Avg95l'],color=c_, alpha=0.75, lw=0.4)
    ax.fill_between(d['mediamovilx'], d['Avg95l'],d['Avg95u'], color=c_, alpha=0.15)

def distanceCI(d, arg, ax, chrom):
    max_x = arg['contigs'][chrom]
    z95=abs(st.norm.ppf(.025/arg['n_markers']))
    d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
    d['DP1']=d['DPdom_1']+d['DPrec_1']
    d['DP2']=d['DPdom_2']+d['DPrec_2']
    d['prodx']=d['POS'] * d['DP']

    d['CI95u']=d['DELTA'] + z95 * np.sqrt(d['SNPidx2']*(1-d['SNPidx2'])/d['DP2'] + d['SNPidx1']*(1-d['SNPidx1'])/d['DP1'])
    d['CI95l']=d['DELTA'] - z95 *np.sqrt(d['SNPidx2']*(1-d['SNPidx2'])/d['DP2'] + d['SNPidx1']*(1-d['SNPidx1'])/d['DP1'])

    d['pCIU95']=d['CI95u'] * d['DP']
    d['pCIL95']=d['CI95l'] * d['DP']

    interval_ranges = [x for x in range(0, max_x, arg['--distance-avg'])]
    if interval_ranges[-1] < max_x:
        interval_ranges.append(max_x)
    interval_labels = [f'{start}-{end-1}' for start, end in zip(interval_ranges, interval_ranges[1:])]
    d['Interval'] = pd.cut(d['POS'], bins=interval_ranges, labels=interval_labels)

    res = d.groupby('Interval').agg({'prodx': sum, 'pCIU95': sum, 'pCIL95':sum, 'DP': sum})
    res['POSx'] = (res['prodx']/res['DP'])
    res['Avg95u'] = (res['pCIU95']/res['DP'])
    res['Avg95l'] = (res['pCIL95']/res['DP'])
    res = res.dropna()
    c_ = arg['--palette']['ci']
    ax.plot(res['POSx'], res['Avg95u'],color=c_, alpha=0.75, lw=0.4)
    ax.plot(res['POSx'], res['Avg95l'],color=c_, alpha=0.75, lw=0.4)
    ax.fill_between(res['POSx'], res['Avg95l'],res['Avg95u'], color=c_, alpha=0.15)

def plot_avg(d, arg, ax, field):
    if field == 'SNPidx2':
        d['DP']=d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['SNPidx2']
    elif field == 'SNPidx1':
        d['DP']=d['DPdom_1']+d['DPrec_1']
        c_=arg['--palette']['SNPidx1']
    elif field == 'MAX_SNPidx2':
        if 'SNPidx2' not in arg['--fields']:
            d['DP']=d['DPdom_1']+d['DPrec_1']
        else:
            d['DP']=d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['MAX_SNPidx2']
    elif field == 'log10PVALUE':
        d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['log10PVALUE']
    elif field == 'ED100_4':
        d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['log10PVALUE']
    elif field == 'DELTA':
        d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['DELTA']
    elif field == 'G':
        d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['mvg']
    elif field == 'delta2':
        d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['SNPidx1']
    d['prody']=d[field] * d['DP']
    d['avgy']=(d['prody'].rolling(arg['--moving-avg']).sum() / \
               d['DP'].rolling(arg['--moving-avg']).sum())
    d['prodx']=d['POS'] * d['DP']
    d['avgx']=(d['prodx'].rolling(arg['--moving-avg']).sum() / \
               d['DP'].rolling(arg['--moving-avg']).sum())
    
    ax.plot(d['avgx'], d['avgy'], c=c_, lw=2)#lw=0.75

def plotDistanceAvg(d, chrom, arg, ax, field):
    if field == 'SNPidx2':
        d['DP']=d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['SNPidx2']
    elif field == 'SNPidx1':
        d['DP']=d['DPdom_1']+d['DPrec_1']
        c_=arg['--palette']['SNPidx1']
    elif field == 'MAX_SNPidx2':
        if 'SNPidx2' not in arg['--fields']:
            d['DP']=d['DPdom_1']+d['DPrec_1']
        else:
            d['DP']=d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['MAX_SNPidx2']
    elif field == 'log10PVALUE':
        d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['log10PVALUE']
    elif field == 'ED100_4':
        d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['log10PVALUE']
    elif field == 'DELTA':
        d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['DELTA']
    elif field == 'G':
        d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['mvg']
    elif field == 'delta2':
        d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
        c_=arg['--palette']['SNPidx1']
    max_x = arg['contigs'][chrom]
    d['prodx'] = d['DP']*d['POS']
    d['prody'] = d['DP']*d[field]
    interval_ranges = [x for x in range(0, max_x, arg['--distance-avg'])]
    if interval_ranges[-1] < max_x:
        interval_ranges.append(max_x)
    interval_labels = [f'{start}-{end-1}' for start, end in zip(interval_ranges, interval_ranges[1:])]
    d['Interval'] = pd.cut(d['POS'], bins=interval_ranges, labels=interval_labels)

    res = d.groupby('Interval').agg({'prodx': sum, 'prody': sum, 'DP': sum})
    res['POSx'] = (res['prodx']/res['DP'])
    res['VALy'] = (res['prody']/res['DP'])
    res = res.dropna()

    ax.plot(res['POSx'], res['VALy'], c=c_, lw=2)


def plot_avg_qtl_SNPidx(d, arg, ax):
    color = arg['--palette']['mvg']
    d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
    d['DP1']=d['DPdom_1']+d['DPrec_1']
    d['DP2']=d['DPdom_2']+d['DPrec_2']
    d['productx']=d['POS'] * d['DP']
    d['mediamovilx']=(d['productx'].rolling(arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())
    #SNPidx1
    d['producty1']=d['SNPidx1'] * d['DP1']
    d['mediamovily1']=(d['producty1'].rolling(arg['--moving-avg']).sum() / d['DP1'].rolling(arg['--moving-avg']).sum())
    ax[0].plot(d['mediamovilx'], d['mediamovily1'], c=arg['--palette']['SNPidx1'], lw=2)
    #SNPidx2
    d['producty2']=d['SNPidx2'] * d['DP2']
    d['mediamovily2']=(d['producty2'].rolling(arg['--moving-avg']).sum() / d['DP2'].rolling(arg['--moving-avg']).sum())
    ax[1].plot(d['mediamovilx'], d['mediamovily2'], c=arg['--palette']['SNPidx2'], lw=2)
    # D-SNPidx
    d['productyD']=d['DELTA'] * d['DP']
    d['mediamovilD']=(d['productyD'].rolling(arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())
    ax[2].plot(d['mediamovilx'], d['mediamovilD'], c=arg['--palette']['DELTA'], lw=2)


def distanceSNPidx(d, arg, ax, chrom):
    max_x = arg['contigs'][chrom]
    d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
    d['DP1']=d['DPdom_1']+d['DPrec_1']
    d['DP2']=d['DPdom_2']+d['DPrec_2']
    d['prodx']=d['POS'] * d['DP']
    #SNPidx1
    d['prody1']=d['SNPidx1'] * d['DP1']
    #SNPidx2
    d['prody2']=d['SNPidx2'] * d['DP2']
    # D-SNPidx
    d['prodyD']=d['DELTA'] * d['DP']
    interval_ranges = [x for x in range(0, max_x, arg['--distance-avg'])]
    if interval_ranges[-1] < max_x:
        interval_ranges.append(max_x)
    interval_labels = [f'{start}-{end-1}' for start, end in zip(interval_ranges, interval_ranges[1:])]
    d['Interval'] = pd.cut(d['POS'], bins=interval_ranges, labels=interval_labels)
    res = d.groupby('Interval').agg({'prodx': sum, 'prody1': sum, 'prody2': sum, 'prodyD': sum,'DP1': sum, 'DP2': sum, 'DP':sum})

    res['POSx'] = (res['prodx']/res['DP'])
    res['avg1'] = (res['prody1']/res['DP1'])
    res['avg2'] = (res['prody2']/res['DP2'])
    res['avgD'] = (res['prodyD']/res['DP'])
    res = res.dropna()
    ax[0].plot(res['POSx'], res['avg1'], c=arg['--palette']['SNPidx1'], lw=2)
    ax[1].plot(res['POSx'], res['avg2'], c=arg['--palette']['SNPidx2'], lw=2)
    ax[2].plot(res['POSx'], res['avgD'], c=arg['--palette']['DELTA'], lw=2)

def plot_boost(d, arg, ax, max_x):
    # Mediamovil Boost
    if 'SNPidx2' not in arg['--fields']:
        d['DP']=d['DPdom_1']+d['DPrec_1']
    else:
        d['DP']=d['DPdom_2']+d['DPrec_2']
    d['BOOST']=d['BOOST'] * arg['lim']
    d['prodboost']=d['BOOST'] * d['DP']
    d['medboost']=(d['prodboost'].rolling(arg['--boost']).sum() / \
                   d['DP'].rolling(arg['--boost']).sum())
    d['prodboostx']=d['POS'] * d['DP']
    d['medboostx']=(d['prodboostx'].rolling(
        arg['--boost']).sum() / d['DP'].rolling(arg['--boost']).sum())
    ax2=ax.twinx()
    ax2.spines['top'].set_visible(False)
    c_=arg['--palette']['BOOST']
    ax2.plot(d['medboostx'], d['medboost'], c=c_, lw=1.25, linestyle='dashed')
    ax2.set(xlim=(0, max_x), ylim=(0, 1))
    ax2.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
    ax2.set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=8)
    ax2.set_ylabel(ylabel='Boost', fontsize=12, rotation=90, labelpad=15)

def plotDistanceBoost(d, arg, ax, max_x, chrom):
    # Mediamovil Boost
    if 'SNPidx2' not in arg['--fields']:
        d['DP']=d['DPdom_1']+d['DPrec_1']
    else:
        d['DP']=d['DPdom_2']+d['DPrec_2']
    d['BOOST']=d['BOOST'] * arg['lim']
    max_x = arg['contigs'][chrom]
    d['prodx'] = d['DP']*d['POS']
    d['prody'] = d['DP']*d['BOOST']
    interval_ranges = [x for x in range(0, max_x, arg['--distance-avg'])]
    if interval_ranges[-1] < max_x:
        interval_ranges.append(max_x)
    interval_labels = [f'{start}-{end-1}' for start, end in zip(interval_ranges, interval_ranges[1:])]
    d['Interval'] = pd.cut(d['POS'], bins=interval_ranges, labels=interval_labels)
    res = d.groupby('Interval').agg({'prodx': sum, 'prody': sum, 'DP': sum})
    res['POSx'] = (res['prodx']/res['DP'])
    res['VALy'] = (res['prody']/res['DP'])
    res = res.dropna()
    ax2=ax.twinx()
    ax2.spines['top'].set_visible(False)
    c_=arg['--palette']['BOOST']
    ax2.plot(res['POSx'], res['VALy'], c=c_, lw=1.25, linestyle='dashed')
    ax2.set(xlim=(0, max_x), ylim=(0, 1))
    ax2.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
    ax2.set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=8)
    ax2.set_ylabel(ylabel='Boost', fontsize=12, rotation=90, labelpad=15)

def plot_ED100_4(d, arg, ax, max_x, max_y):
    ax2 = ax.twinx()
    c_ = arg['--palette']['mvg']
    ax2.set(xlim=(0, max_x), ylim=(0, max_y))
    ax2.tick_params(axis='y', which='major', labelsize=8)
    ax2.set_ylabel(ylabel='ED $\mathregular{100^4}$  $\mathregular{x10^8}$ ', fontsize=10, rotation=90, labelpad=15)
    ax2.plot(d['POS'], d['ED100_4'], c=c_, lw=2)
    ax2.spines['top'].set_visible(False)

def create_caption(arg, res_tit):
    f=open(arg['captions_dir'] + res_tit + '.txt', 'a')
    return f

def write_caption(f, text, arg):
    for sentence in text:
        f.write(sentence+' ')
    if '--window' in arg.keys():
        f.write('The dots correspond to genomic regions defined by non-overlapping bins of {} consecutive markers'.format(str(arg['--window'])))
    if '--distance' in arg.keys():
        f.write('The dots correspond to bins of all markers within a range of {} pb'.format(str(arg['--distance'])))
    f.write('\n'*2)
    f.close()

def Delta2_Vertical_graph(df, arg):
    typ=arg['--output-type']
    max_x=arg['max_lenght']
    
    ticks_y=[-1,-0.5,0,0.5, 1]
    lim_y=(-1, 1)
    #t1,rt1,_=arg['titles'][7]
    #t2,rt2,_=arg['titles'][8]
    labs_list = list()
    f1_name = list()
    f2_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig, ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            d['delta2'] = d['DPdom_1']/(d['DPdom_1']+d['DPrec_1']) - d['DPdom_2']/(d['DPdom_2']+d['DPrec_2'])
            y=d[['delta2']]
            #AF1&2
            ax[i].scatter(x, y, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
            ax[i].set(xlim=(0, max_x), ylim=lim_y)
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='Allele Frequency', fontsize=12, rotation=90, labelpad=15)
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=16, rotation=0, x = -0.12, y=0.85) 
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))  
            ax[i].set_yticks(ticks=ticks_y)
            ax[i].set_yticklabels(labels=ticks_y, fontsize=8)
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], 'delta2')
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            if len(chrom) == 1:
                ax[1].remove()
            if i == len(chrom)-1 or len(chrom) == 1:
                ax[i].tick_params(labelbottom=True)
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=8, color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)

        fig.subplots_adjust(hspace=0.1)
        filename='delta2multiV'+typ
        filename=check_save(arg, filename)
        f1_name.append(filename)
        plt.savefig(arg['--outdir']+filename, dpi=arg['DPI'])
        plt.close()