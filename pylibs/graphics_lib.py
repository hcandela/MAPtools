from pylibs.analysis_lib import *
from pylibs.constants import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
import scipy.stats as st
import sys
import os
import io
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

def read_header_plot(arg):
    arg['--contigs'] = dict()
    with open(arg['--input'],'r') as handle:
        for line in handle:
            if line.startswith('##maptools_mbsCommand='):
                arg['type'] = 'mbs'

                names = ['-r', '--ref-genotype','-d','--data']
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
                        arg['--data'] = data
                if 'D' in data and 'R' in data:
                    arg['pool1'] = 'dominant'
                    arg['pool2'] = 'recessive'
                else:
                    arg['pool2'] = 'recessive'

                arg['--ref-genotype'] = True if True in result else False
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
                arg['type'] = 'qtl'
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
                        arg['--data'] = data
                arg['pool1'] = 'high'
                arg['pool2'] = 'low'
                arg['--ref-genotype'] = True if True in result else False

            if line.startswith('##contig=<'):
                s = line[line.find('<')+1:line.rfind('>')].split(',')
                id = s[0].split('=')[1]
                length = s[1].split('=')[1]
                arg['--contigs'][id] = int(length)

            if line.startswith('#CHROM'):
                arg['header'] = line.rstrip().split('\t')
                if arg['type'] == 'qtl':
                    translate = {'HIGH':'DOM','LOW':'REC','DPhigh_1':'DPdom_1','DPlow_1':'DPrec_1','DPhigh_2':'DPdom_2','DPlow_2':'DPrec_2'}
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
    if not arg['--input']:
        print(__doc__, end='\n', file=sys.stderr)
        sys.exit()
    arg['--input'] = os.path.abspath(arg['--input'])
    inp_f = arg['--input']
    if not inp_f:
        print(__doc__, end='\n', file=sys.stderr)
        sys.exit()
    try:
        f = open(os.path.abspath(arg['--input']), 'r')
    except FileNotFoundError:
        print('Error: The input file {} does not exist'.format(inp_f), file=sys.stderr)
        sys.exit()
    wd = os.getcwd()
    try:
        absolut_path = os.path.abspath(arg['--outdir'])
        os.makedirs(absolut_path)
    except FileExistsError:
        print('Warning: the output directory already exists', file=sys.stderr)
        pass

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
        print('Error: Options -A and -W are incompatible.', file=sys.stderr)
        sys.exit()
    return arg

def load_dataframe_plotting(arg):
    read_settings(arg)

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
        df1 = df.copy()
        df1['#CHROM'] = pd.Categorical(df1['#CHROM'], categories=chroms_, ordered=True)
        df = df1.sort_values(by=['#CHROM','POS'])
        df = df.reset_index()

    if arg['alias'] != False:
        if set(arg['alias']).issubset(arg['--chromosomes']):
            arg['--contigs'] = {arg['alias'].get(k,k): v for k, v in arg['--contigs'].items()}
            arg['--chromosomes'] = [arg['alias'].get(c, c) for c in arg['--chromosomes']]
            df['#CHROM'] = df['#CHROM'].map(arg['alias'])
        else:
            print('Warning: The chromosome names do not match those in the alias list. Names may not be substituted.', file=sys.stderr)

    if len(arg['--chromosomes']) > 1 and arg['--multi-chrom']:
        arg['--multi-chrom'] = True
    else:
        arg['--multi-chrom'] = False
    arg = check_chroms(arg)
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
    try:
        float(arg['--dot-size'])
    except ValueError:
        print('Warning: The --dot-size argument cannot be a literal.', file=sys.stderr)
        arg['--dot-size'] = 1
    
    if float(arg['--dot-size']) < 0:
        print('Warning: the dot size cannot be less than zero.', file=sys.stderr)   
        arg['DOT_SIZE'] = 0
    else:
        arg['DOT_SIZE'] = float(arg['--dot-size'])
    

    try:
        float(arg['--line-thickness'])
    except ValueError:
        print('Warning: The --line-thickness argument cannot be a literal.', file=sys.stderr)
        arg['--line-thickness'] = 2

    if float(arg['--line-thickness']) < 0:
        print('Warning: the line width cannot be less than zero.', file=sys.stderr)   
        arg['LINE_W'] = 1
    else:
        arg['LINE_W'] = float(arg['--line-thickness'])

    try:
        int(arg['--DPI'])
    except ValueError:
        print('Warning: The --DPI argument cannot be a literal or float.', file=sys.stderr)
        arg['--DPI'] = 600
    
    if int(arg['--DPI']) < 0:
        print('Warning: dots per inch argument cannot be less than zero.', file=sys.stderr)
    else:
        arg['DPI'] = int(arg['--DPI'])

    arg = checkPlottingOptions(arg,df)
    return arg, df

def checkPlottingOptions(arg,df):
    arg['titles'] = TITLES
    arg['lines'] = LINES
    arg['ED100'] = False
    if arg['type'] == 'mbs':
        if arg['--boost'] == True:
            if arg['--moving-avg'] == False and arg['--distance-avg'] == False:
                arg['--boost'] = False
                print('Warning: Activate the -A or -W options to use boost.', file=sys.stderr)

    if arg['--all']:
        if 'log10PVALUE' in arg['--fields']:
            arg['--pvalue'] = True
        else:
            arg['--pvalue'] = False
        if 'SNPidx1' in arg['--fields']:
            arg['--allele-freq-1'] = True
        if 'SNPidx2' in arg['--fields']:
            arg['--allele-freq-2'] = True
        # Activaba Combine cuando usabas -a en todos los gráficos
        if 'SNPidx1' in arg['--fields'] and 'SNPidx2' in arg['--fields'] and arg['--combine'] == True:
            arg['--combine'] = True
        else:
            arg['--combine'] = False
        if len(arg['--chromosomes']) > 1:
            arg['--multi-chrom'] = True
        if 'MAX_SNPidx2' in arg['--fields'] and arg['type'] == 'mbs':
            arg['--max-allele-freq2'] = True
        else:
            arg['--max-allele-freq2'] = False

        if 'DELTA' in arg['--fields']:
            arg['--delta'] = True
            #arg['--ci95'] = True
        if 'ED' in arg['--fields']:
            arg['--euclidean-distance'] = True
            ed100 = get_ED100_4(df, arg, RANG)
            res = list()
            for ch in arg['--chromosomes']:
                d = df[df['#CHROM'] == ch]
                if len(d) < RANG:
                    res.append(False)
                    print(f'Warning: there are not enough markers in {ch} to calculate ED100', file=sys.stderr)
                else:
                    res.append(True)
            if True in res:
                arg['ED100'] = ed100
            else:
                arg['ED100'] = False
        if 'G' in arg['--fields']:
            arg['--g-statistic'] = True
        if arg['--pvalue'] and arg['--delta'] and arg['--euclidean-distance'] and arg['--g-statistic'] and arg['--allele-freq-1'] and arg['--allele-freq-2']:
            arg['--comb-statistics'] = True
    else:
        if arg['--combine']:
            if 'SNPidx1' in arg['--fields'] and 'SNPidx2' in arg['--fields'] and (arg['--moving-avg'] != False or arg['--distance-avg'] != False):
                arg['--combine'] = True
            else:
                print('Warning: is not possible make combined graphics with your input data', file=sys.stderr)
                arg['--combine'] = False
        
        if arg['--pvalue'] == True and 'log10PVALUE' not in arg['--fields']:
            print('Warning: is not possible make P-VALUE graphics with your input data', file=sys.stderr)
            arg['--pvalue'] = False
        if arg['--allele-freq-1'] == True and 'SNPidx1' not in arg['--fields']:
            print('Warning: is not possible make phased Allele Frequency graphics with your input data. Please use -M or -a option', file=sys.stderr)
            arg['--allele-freq-1'] = False
        if arg['--allele-freq-2'] == True and 'SNPidx2' not in arg['--fields']:
            print('Warning: is not possible make Allele Frequency graphics for this pool. Please use -R, -M or -a option', file=sys.stderr)
            arg['--allele-freq-2'] = False
        if arg['--max-allele-freq2'] == True and arg['type'] == 'mbs' and 'MAX_SNPidx2' not in arg['--fields']:
            print('Warning: is not possible make phased MAX Allele Frequency graphics with your input data.', file=sys.stderr)
            arg['--max-allele-freq2'] = False
        if arg['--max-allele-freq2'] == True and arg['type'] == 'qtl':
            print('Warning: is not possible make MAX Allele Frequency graphics with your QTL-Seq input data.', file=sys.stderr)
            arg['--max-allele-freq2'] = False
        if arg['--delta'] == True and 'DELTA' not in arg['--fields']:
            print('Warning: is not possible make phased SNP-idx graphics with your input data. Please use -a.', file=sys.stderr)
            arg['--delta'] = False
        if arg['--ci95'] == True and 'DELTA' not in arg['--fields']:
            print('Warning: is not possible to calculate confidense interval without DELTA field. Please use -a.', file=sys.stderr)
            arg['--ci95'] = False
        if arg['--euclidean-distance'] == True:
            if 'ED' in arg['--fields']:
                ed100 = get_ED100_4(df, arg, RANG)
                res = list()
                for ch in arg['--chromosomes']:
                    d = df[df['#CHROM'] == ch]
                    if len(d) < RANG:
                        res.append(False)
                        print(f'Warning: there is not enough markers in {ch} to calculate ED100', file=sys.stderr)
                    else:
                        res.append(True)
                if True in res:
                    arg['ED100'] = ed100
                else:
                    arg['ED100'] = False
            else:
                print('Warning: is not possible make Euclidean Distance graphics with your input data', file=sys.stderr)
                arg['--euclidean-distance'] = False
        if arg['--g-statistic'] == True and 'G' not in arg['--fields']:
            print('Warning: is not possible make G-statistic graphics with your input data', file=sys.stderr)
            arg['--g-statistic'] = False
        if arg['--comb-statistics'] == True:
            if all(x in arg['--fields'] for x in['SNPidx1', 'SNPidx2', 'DELTA', 'ED', 'G', 'log10PVALUE']):
                arg['--comb-statistics'] = True
            else:
                arg['--comb-statistics'] = False

            if not isinstance(arg['ED100'],pd.DataFrame):
                ed100 = get_ED100_4(df, arg, RANG)
                res = list()
                for ch in arg['--chromosomes']:
                    d = df[df['#CHROM'] == ch]
                    if len(d) < RANG:
                        res.append(False)
                        print(f'Warning: there is not enough markers in {ch} to calculate ED100', file=sys.stderr)
                    else:
                        res.append(True)
                if True in res:
                    arg['ED100'] = ed100
                else:
                    arg['ED100'] = False

    if arg['type'] == 'mbs':
        if arg['--all'] == False and arg['--max-allele-freq2'] and arg['--pvalue'] == False and arg['--allele-freq-1'] == False and arg['--allele-freq-2'] and arg['--delta'] == False and arg['--comb-statistics'] == False and arg['--g-statistic'] == False and arg['--euclidean-distance'] == False:
            print('Warning: any graphic was selected. Please use -a option to make all possible graphics with your data.', file=sys.stderr)
    else:
        if arg['--all'] == False and arg['--pvalue'] == False and arg['--allele-freq-1'] == False and arg['--allele-freq-2'] and arg['--delta'] == False and arg['--comb-statistics'] == False and arg['--g-statistic'] == False and arg['--euclidean-distance'] == False:
            print('Warning: any graphic was selected. Please use -a option to make all possible graphics with your data.', file=sys.stderr)
    return arg

def load_dataframe_merge(arg):
    inp_f = arg['--input']
    inp_ext = inp_f.split('.')[-1]
    if inp_ext == 'csv':
        sep_ = ','
    else:
        sep_ = '\t'
    df = pd.read_csv(inp_f, sep=sep_, dtype=fields, names=arg['header2'], comment='#')
    arg['--fields'] = list(df.columns)
    arg['header']=arg['--fields']
    if arg['--chromosomes'] == 'all':
        arg['--chromosomes'] = list(df['#CHROM'].unique())
    else:
        arg['--chromosomes'] = arg['--chromosomes'].split(',')
    arg = check_chroms(arg)
    return arg, df
    
def check_merge(arg,__doc__):
    arg['head'] = dict()
    arg['head'][8] = list()
    arg['head'][9] = list()
    global fields
    inp_f = arg['--input']
    if not inp_f:
        print(__doc__,end='\n', file=sys.stderr)
        sys.exit()
    arg['--input'] = os.path.abspath(arg['--input'])
    inp_f = os.path.abspath(arg['--input'])
    try:
        f = open(arg['--input'], 'r')
    except FileNotFoundError:
        print('Error: The input file {} does not exist'.format(inp_f), file=sys.stderr)
        sys.exit()

    arg['spacer'] = '\t'
    arg['--output-type'] = '.txt'

    if arg['--output'] != None:
        outdir = os.path.dirname(arg['--output']) or os.getcwd()
        absolut_path = os.path.abspath(outdir)
        try:
            os.makedirs(absolut_path)
        except FileExistsError:
            print('Warning: the output directory already exists', file=sys.stderr)
            pass
        arg['filename'] = arg['--output'].split('/')[-1]
        arg['outdir'] = absolut_path
    
        arg['--output'] = arg['outdir'] + '/'+ arg['filename']
        if arg['--output'] != None:
            arg['--output'] = check_save_an(arg, arg['filename'])
        else:
            arg['--output'] = None

    arg['lim'] = 1e-90

    if arg['--window'] != None and arg['--distance'] != None:
        print('Error: Options -w and -D are incompatible.', file=sys.stderr)
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


def read_settings(arg):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    json_file_path = os.path.join(script_dir,'..','settings.json')
    try:
        with open(json_file_path) as json_file:
            data = json.load(json_file)
    except FileNotFoundError:
        print('Error: Please put the file \"settings.json\" in the MAPtools folder.', file=sys.stderr)
        sys.exit()

    s = data['settings']
    for key,dicc in s.items():
        for k,v in dicc.items():
            if v == 'None':
                dicc[k] = None
    arg['sets'] = s
    
    a = data['alias']
    if a == dict():
        arg['alias'] = False
    else:
        arg['alias'] = a
    
    p = data['palettes']
    palette = arg['--palette']
    if palette not in p.keys():
        print(f'Warning: {palette} is not one of the available palette names.\nSwitching to standard', file=sys.stderr)
        palette = 'standard'
    
    for key,val in p[palette].items():
        if val == list():
            palette = 'standard'
            break
    arg['--palette'] = {k: v[0] for k, v in p[palette].items()}
    arg['color_names'] = {k: v[1] for k, v in p[palette].items()}


def check_save(arg, file_name):
    typ = arg['--output-type']
    if os.path.isfile(arg['sub_folder']+file_name):

        expand = 0
        while True:
            expand += 1
            nw_file_name = file_name.split(typ)[0] + '_' + str(expand) + typ
            if os.path.isfile(arg['sub_folder']+nw_file_name):
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


def grouped_by(df, arg):
    w = arg['--window']
    chrom = arg['--chromosomes']
    spacer = arg['spacer']
    fsal = arg['fsal']
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
            if 'DPdom_1' not in arg['--fields'] and 'DPrec_1' not in arg['--fields']:
                rAt, rBt, rPost, rTt = 0, 0, 0, 0
            else:
                rAt, rBt, rCt, rDt, rPost, rTt = 0, 0, 0, 0, 0, 0
            for i in g:  # for each index in group
                if 'DPdom_1' not in arg['--fields'] and 'DPrec_1' not in arg['--fields']:
                    rA = d.loc[i].DPdom_2
                    rB = d.loc[i].DPrec_2
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
            if 'DPdom_1' in arg['--fields'] and 'DPrec_1' in arg['--fields']:
                res += [int(rCt), int(rDt)]
                if 'SNPidx1' in arg['--fields'] and 'SNPidx2' in arg['--fields']:
                    rSNPidx1 = rBt/(rAt+rBt)
                    rSNPidx2 = rDt/(rCt+rDt)
                    res += [rSNPidx1, rSNPidx2]
                if 'MAX_SNPidx2' in arg['--fields']:
                    rMAX_SNPidx2 = (max(rDt, rCt))/(rDt+rCt)
                    rSNPidx1 = rBt/(rAt+rBt)
                    rSNPidx2 = rDt/(rCt+rDt)
                    fisher = LogFisher(rAt, rBt, rCt, rDt)
                    res += [rMAX_SNPidx2]
                if 'FISHER' in arg['--fields']:
                    fisher = LogFisher(rAt, rBt, rCt, rDt)
                    res += [fisher]

                pval, log10pval = pvalor2(rAt, rBt, rCt, rDt)
                rDELTA = rSNPidx2 - rSNPidx1
                v1 = (rCt/(rCt + rDt), rDt/(rCt + rDt))
                v2 = (rAt/(rAt + rBt), rBt/(rAt + rBt))
                ed = distance.euclidean(v1,v2)
                g = Gstatic(rAt, rBt, rCt, rDt)
                res += [pval, log10pval, rDELTA, ed, g]
            else:
                if 'MAX_SNPidx2' in arg['--fields']:
                    rMAX_SNPidx2 = (max(rAt, rBt))/(rAt+rBt)
                    res += [rMAX_SNPidx2]
                elif 'SNPidx2' in arg['--fields']:
                    rSNPidx2 = rBt/(rAt+rBt)
                    res += [rSNPidx2]

            n_line = spacer.join(str(field) for field in res) + '\n'
            write_line(n_line, fsal)




def get_ED100_4(df, arg, rang):
    #chrom = [ch for ch in arg['contigs'].keys()]
    chrom = arg['--chromosomes']
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
    ed100 = df.groupby('#CHROM').filter(lambda x: len(x) > 100)
    ed100 = ed100[['#CHROM','POS', 'ED','ED100_4']]
    return ed100


def EDmonoPlot(df, arg):
    sets = arg['sets']['mono']
    chrom = arg['--chromosomes']
    typ = arg['--output-type']
    ed100 = arg['ED100']
    if isinstance(ed100,pd.DataFrame):
        max_y = max(ed100['ED100_4'])
    t,rt=arg['titles'][13]
    for ch in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[ch]]
        d.index = np.arange(len(d))
        max_x = int(arg['contigs'][chrom[ch]])
        x = d[['POS']]
        y = d[['ED']]
        fig, ax = plt.subplots(figsize=(sets['figwidth'], sets['figheight']))
        ax.scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax.set(xlim=(0, max_x), ylim=(0, 1.5))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=sets['xticksSIZE'])
        ax.tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Genomic position (pb)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax.set_ylabel(ylabel='EDm', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
        ax.set_yticks(ticks=[0, 0.5, 1, 1.5])
        ax.set_yticklabels(labels=[0, 0.5, 1, 1.5], fontsize=sets['yticksSIZE'])
        ax.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax.spines['top'].set_visible(False)
        for spine in ax.spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if isinstance(ed100,pd.DataFrame):
            plot_ED100_4(arg, ax, max_x, max_y, ed100, chrom[ch], sets)
            
        rtch = rt.format(chrom[ch])
        filename = rtch + typ
        filename = check_save(arg, filename)
        plt.tight_layout()
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
        del d
        cap = list()
        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t.format(chrom[ch]))
            cap.append(arg['lines'][23].format(arg['color_names']['dots']))
            if isinstance(ed100,pd.DataFrame):
                cap.append(arg['lines'][24].format(arg['color_names']['mvg']))
            write_caption(f,cap,arg)

def GmonoPlot(df, arg):
    sets = arg['sets']['mono']
    chrom = arg['--chromosomes']
    typ = arg['--output-type']
    min_y, max_y =min(df['G']), max(df['G'])
    t,rt = arg['titles'][15]
    for ch in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[ch]]
        max_x = int(arg['contigs'][chrom[ch]])
        x = d[['POS']]
        y = d[['G']]
        fig, ax = plt.subplots(figsize=(sets['figwidth'], sets['figheight']))
        ax.scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax.set(xlim=(0, max_x), ylim=(min_y, max_y))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=sets['xticksSIZE'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Chromosomal position (pb)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax.set_ylabel(ylabel='G-statistic', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
        ax.tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax.spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax, 'G')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[ch], arg, ax, 'G')
        rtch = rt.format(chrom[ch])
        filename = rtch + typ
        filename = check_save(arg, filename)
        plt.tight_layout()
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
        cap = list()

        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t.format(chrom[ch]))
            cap.append(arg['lines'][25].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][27].format(arg['color_names']['mvg'], str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][26].format(arg['color_names']['mvg'], str(arg['--distance-avg'])))
            write_caption(f,cap,arg)
        del d


def AF_manhattan_plot(df, arg, g_type):
    sets = arg['sets']['manhattan']
    typ = arg['--output-type']
    #labs_list = list()
    f_name = list()
    ylab = 'SNP-index'
    boost = False
    multi = (False, True)
    ticks_y=[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    arg, ticks_y, lim_y, ylab, t, rt, key2, key = chooseOptionsAFPlots(arg, multi, g_type, ticks_y, lim_y, ylab)
    chrom = arg['--chromosomes']
    if len(arg['--chromosomes'])*sets['chromMinWIDTH'] <= sets['maxWIDTH']:
        af = len(arg['--chromosomes'])*sets['chromMinWIDTH']
    else:
        af = sets['maxWIDTH']
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, sets['figheight']))

    for i in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[i]]
        max_x = int(arg['contigs'][chrom[i]])
        x = d[['POS']]
        y = d[[g_type]]
        ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[i].set(xlim=(0, max_x), ylim=lim_y)
        ax[i].set_xticks([])
        ax[i].set_xlabel(xlabel=f'{chrom[i]}', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        #labs_list.append('Chromosome {}'.format(chrom[i]))
        ax[i].set_ylabel(ylabel=ylab,fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'])
        ax[i].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[i].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[i].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[i], g_type)
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[i], g_type)

        if arg['type'] == 'mbs':
            if arg['--boost'] == True and arg['--moving-avg'] != False:
                keyB = 36
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plot_boostManhattan(d, arg, ax[i], max_x, chrom, i, sets, ticks_y, g_type)
                    boost = True
            if arg['--boost'] == True and arg['--distance-avg'] != False:
                keyB = 35
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plotDistanceBoostManhattan(d, arg, ax[i], max_x, chrom, i, sets, ticks_y, g_type)
                    boost = True
        if arg['--ci95'] and g_type == 'DELTA':
            if arg['--moving-avg'] != False:
                calc_ci(d, arg, ax[i])
            if arg['--distance-avg'] != False:
                distanceCI(d, arg, ax[i], chrom[i])
            ax[i].axhline(y=0, color='black', linestyle='dashed', linewidth=0.75)

        ax[0].set_yticks(ticks=ticks_y)
        ax[0].set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        if len(chrom) == 1:
            ax[1].remove()
        if i > 0:
            ax[i].axes.get_yaxis().set_visible(False)
        del d
    #fig.tight_layout()
    fig.subplots_adjust(wspace=sets['wspace'],bottom=sets['bottom'])
    filename = rt + typ
    filename = check_save(arg, filename)
    f_name.append(filename)
    plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
    plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        cap.append(arg['lines'][key2].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][key].format(arg['color_names'][g_type], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][key].format(arg['color_names'][g_type], str(arg['--distance-avg'])))
        if arg['type'] == 'mbs':
            if boost == True:
                if arg['--moving-avg'] != False:
                    cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--moving-avg'])))
                if arg['--distance-avg'] != False:
                    cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--distance-avg'])))
        if arg['--ci95'] and g_type == 'DELTA' and (arg['--moving-avg'] != False or arg['--distance-avg'] != False):
                cap.append(arg['lines'][22].format(str(arg['n_markers'])))
        write_caption(f, cap,arg)


def AFCombinedManhattanPlot(df, arg):
    sets = arg['sets']['manhattan']
    typ = arg['--output-type']
    #labs_list = list()
    f_name = list()
    ylab = 'SNP-index'
    ticks_y=[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    t,rt = arg['titles'][4]

    chrom = arg['--chromosomes']
    if len(arg['--chromosomes'])*sets['chromMinWIDTH'] <= sets['maxWIDTH']:
        af = len(arg['--chromosomes'])*sets['chromMinWIDTH']
    else:
        af = sets['maxWIDTH']
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, sets['figheight']))

    for i in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[i]]
        max_x = int(arg['contigs'][chrom[i]])
        x = d[['POS']]
        #y = d[[g_type]]
        #ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'],  clip_on=False clip_on=False)
        ax[i].set(xlim=(0, max_x), ylim=lim_y)
        ax[i].set_xticks([])
        ax[i].set_xlabel(xlabel=f'{chrom[i]}', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        #labs_list.append('Chromosome {}'.format(chrom[i]))
        ax[i].set_ylabel(ylabel=ylab,fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'])
        ax[i].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[i].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[i].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        #SNPidx1
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[i], 'SNPidx1')
            if arg['type'] == 'mbs':
                l1,l2 = 3,5
            else:
                l1,l2 = 9,11
        if arg['--distance-avg'] != False:
            if arg['type'] == 'mbs':
                l1,l2 = 2,4
            else:
                l1,l2 = 8,10
            plotDistanceAvg(d, chrom[i], arg, ax[i], 'SNPidx1')

        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[i], 'SNPidx2')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[i], 'SNPidx2')
        ax[i].axhline(y=0.5, color='black', linestyle='dashed', linewidth=0.75)

        ax[0].set_yticks(ticks=ticks_y)
        ax[0].set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        if len(chrom) == 1:
            ax[1].remove()
        if i > 0:
            ax[i].axes.get_yaxis().set_visible(False)
        del d
    #fig.tight_layout()
    fig.subplots_adjust(wspace=sets['wspace'], bottom=sets['bottom'])
    filename = rt + typ
    filename = check_save(arg, filename)
    f_name.append(filename)
    plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
    plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][l2].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
            cap.append(arg['lines'][l1].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][l2].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
            cap.append(arg['lines'][l1].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])))
        write_caption(f, cap,arg)


def pval_manhattan_plot(df, arg):
    sets = arg['sets']['manhattan']
    t,_,rt = arg['titles'][12]
    typ = arg['--output-type']
    max_y = max(-df['log10PVALUE'])*1.05
    #labs_list = list()
    f_name = list()
    chrom = arg['--chromosomes']
    if len(arg['--chromosomes'])*sets['chromMinWIDTH'] <= sets['maxWIDTH']:
        af = len(arg['--chromosomes'])*sets['chromMinWIDTH']
    else:
        af = sets['maxWIDTH']
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, sets['figheight']))
    threshold = -log(0.05/len(df.axes[0]))/log(10)
    for i in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[i]]
        max_x = int(arg['contigs'][chrom[i]])
        x = d[['POS']]
        y = -d[['log10PVALUE']]
        ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[i].set(xlim=(0, max_x), ylim=(0, (max(max_y, threshold)+1)//1))
        ax[i].set_xticks([])
        ax[i].set_xlabel(xlabel=f'{chrom[i]}', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        #labs_list.append('Chromosome {}'.format(chrom[i]))
        ax[i].set_ylabel(ylabel='-log'+r'$_{10}$'+'(p-value)',fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'])
        ax[i].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[i].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[i].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if arg['--moving-avg'] != False:
            d['DP']=d['DPdom_1']+d['DPrec_1']+d['DPdom_2']+d['DPrec_2']
            c_=arg['--palette']['log10PVALUE']
            d['prody']=d['log10PVALUE'] * d['DP']
            d['avgy']=(d['prody'].rolling(arg['--moving-avg']).sum() / \
                d['DP'].rolling(arg['--moving-avg']).sum())
            d['prodx']=d['POS'] * d['DP']
            d['avgx']=(d['prodx'].rolling(arg['--moving-avg']).sum() / \
            d['DP'].rolling(arg['--moving-avg']).sum())
            ax[i].plot(d['avgx'], -d['avgy'], c=c_, lw=arg['LINE_W'])#lw=0.75
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
            ax[i].plot(res['POSx'], -res['VALy'], c=c_, lw=arg['LINE_W'])

        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        if arg['--bonferroni'] and not d.empty:
            max_x_ch = (max(d['POS'])/max_x)
            ax[i].axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        if len(chrom) == 1:
            ax[1].remove()
        if i > 0:
            ax[i].axes.get_yaxis().set_visible(False)
        del d
    #fig.tight_layout()
    fig.subplots_adjust(wspace=sets['wspace'],bottom=sets['bottom'])
    filename = rt + typ
    filename = check_save(arg, filename)
    f_name.append(filename)
    plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
    plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t.format(''))
        cap.append(arg['lines'][12].format(arg['color_names']['dots'],'-'))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][14].format(arg['color_names']['log10PVALUE'],'-', str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][13].format(arg['color_names']['log10PVALUE'], '-', str(arg['--distance-avg'])))
        if arg['--bonferroni']:
            cap.append(arg['lines'][15].format(str(arg['n_markers'])))
        write_caption(f, cap,arg) 

def pval_mono_graph(df, arg):
    sets = arg['sets']['mono']
    chrom=arg['--chromosomes']
    typ=arg['--output-type']
    threshold=log(0.05/len(df.axes[0]))/log(10)
    min_y=min(df['log10PVALUE'])*1.05
    
    for i in range(len(chrom)):
        d=df[df['#CHROM'] == chrom[i]]
        max_x= int(arg['contigs'][chrom[i]])
        t, rt = arg['titles'][11]
        cap = list()
        x=d[['POS']]
        y=d[['log10PVALUE']]
        fig, ax=plt.subplots(figsize=(sets['figwidth'], sets['figheight']))
        ax.scatter(x, y, s=arg['DOT_SIZE'], color=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax.set(xlim=(0, max_x), ylim=((min(min_y, threshold)-1)//1, 0))
        ax.tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=sets['xticksSIZE'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Genomic position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax.set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)', fontsize=sets['ylabSIZE'], labelpad=sets['ylabDIST'])
        ax.tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax.spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax, 'log10PVALUE')
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax, 'log10PVALUE')
        if arg['--bonferroni'] and not d.empty:
            max_x_ch=(max(d['POS'])/max_x)
            ax.axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        rtch = rt.format(chrom[i])

        filename=rtch + typ
        filename=check_save(arg, filename)
        plt.tight_layout()
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
        del d

        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t.format(chrom[i]))
            cap.append(arg['lines'][12].format(arg['color_names']['dots'],''))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][14].format(arg['color_names']['log10PVALUE'], '', str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][13].format(arg['color_names']['log10PVALUE'], '', str(arg['--distance-avg'])))
            if arg['--bonferroni']:
                cap.append(arg['lines'][15].format(str(arg['n_markers'])))
            write_caption(f, cap,arg)



def AF1_AF2_mono_graph(df, arg):
    sets = arg['sets']['mono']
    chrom=arg['--chromosomes']
    typ=arg['--output-type']
    if arg['type'] == 'mbs':
        k1, k2 = 0,1
        ti1, ti2 = 0,1
    else:
        k1,k2 = 6,7
        ti1,ti2 = 5,6
    for i in range(len(chrom)):
        d=df[df['#CHROM'] == chrom[i]]
        t1,rt1,_=arg['titles'][ti1]
        t2,rt2,_=arg['titles'][ti2]
        max_x= int(arg['contigs'][chrom[i]])
        x=d[['POS']]
        y1=d[['SNPidx1']]
        y2=d[['SNPidx2']]
        # AF1
        fig, ax=plt.subplots(figsize=(sets['figwidth'], sets['figheight']))
        ax.scatter(x, y1, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax.set(xlim=(0, max_x), ylim=(0, 1))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=sets['xticksSIZE'])
        ax.tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Genomic position (bp)'.format(chrom[i]), fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax.set_ylabel(ylabel='SNP-index1', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
        ax.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
        ax.set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=sets['yticksSIZE'])
        ax.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax.spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax, 'SNPidx2')
            plot_avg(d, arg, ax, 'SNPidx1')
            if arg['type'] == 'mbs':
                l1,l2 = 3,5
            else:
                l1,l2 = 9,11
        if arg['--distance-avg'] != False:
            if arg['type'] == 'mbs':
                l1,l2 = 2,4
            else:
                l1,l2 = 8,10
            plotDistanceAvg(d, chrom[i], arg, ax, 'SNPidx2')
            plotDistanceAvg(d, chrom[i], arg, ax, 'SNPidx1')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        rtch1 = rt1.format(chrom[i])
        filename=rtch1 + typ
        filename=check_save(arg, filename)
        plt.tight_layout()
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
        cap = list()
        if arg['--captions']:
            f = create_caption(arg, rtch1)
            cap.append(filename)
            cap.append(t1.format(chrom[i]))
            cap.append(arg['lines'][k1].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][l1].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])))
                cap.append(arg['lines'][l2].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][l1].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])))
                cap.append(arg['lines'][l2].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
            write_caption(f,cap,arg) 
        # AF2
        fig, ax=plt.subplots(figsize=(sets['figwidth'], sets['figheight']))
        ax.scatter(x, y2, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax.set(xlim=(0, max_x), ylim=(0, 1))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=sets['xticksSIZE'])
        ax.tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Genomic position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax.set_ylabel(ylabel='SNP-index2', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
        ax.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
        ax.set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=sets['yticksSIZE'])
        ax.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax.spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
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
        plt.tight_layout()
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
        cap = list()
        if arg['--captions']:
            f = create_caption(arg, rtch2)
            cap.append(filename)
            cap.append(t2.format(chrom[i]))
            cap.append(arg['lines'][k2].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][l2].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
                cap.append(arg['lines'][l1].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])))
            if arg['--distance-avg'] != False:
                cap.append(arg['lines'][l2].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
                cap.append(arg['lines'][l1].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])))
            write_caption(f,cap,arg)
        del d


def chooseOptionsAFPlots(arg, multi, g_type, ticks_y, lim_y, ylab):
    key = False
    if arg['type'] == 'mbs':
        if g_type == 'SNPidx1':
            ylab = 'SNP-index1'
            key2 = 0
            if arg['--moving-avg'] != False:
                key = 3
            if arg['--distance-avg'] != False:
                key = 2
            if not multi[0] and not multi[1]:
                t,_,rt=arg['titles'][0]
            else:
                if multi[0]:
                    t,_,rt,_=arg['titles'][2]
                if multi[1]:
                    t,_,_,rt=arg['titles'][2]

        if g_type == 'SNPidx2':
            ylab = 'SNP-index2'
            key2 = 1
            if arg['--moving-avg'] != False:
                key = 5
            if arg['--distance-avg'] != False:
                key = 4
            if not multi[0] and not multi[1]:
                t,_,rt=arg['titles'][1]
            else:
                if multi[0]:
                    t,_,rt,_=arg['titles'][3]
                if multi[1]:
                    t,_,_,rt=arg['titles'][3]

        if g_type == 'MAX_SNPidx2':
            key2 = 32
            if arg['--moving-avg'] != False:
                key = 33
            if arg['--distance-avg'] != False:
                key = 34
            ticks_y=[0.5, 0.75, 1]
            lim_y=(0.5, 1)
            ylab = 'Max SNP-idx2'
            if not multi[0] and not multi[1]:
                t,rt=arg['titles'][18]
            else:
                if multi[0]:
                    t,rt,_=arg['titles'][19]
                if multi[1]:
                    t,_,rt=arg['titles'][19]
    else:
        if g_type == 'SNPidx1':
            ylab = 'SNP-index1'
            key2 = 6
            if arg['--moving-avg'] != False:
                key = 9
            if arg['--distance-avg'] != False:
                key = 8
            if not multi[0] and not multi[1]:
                t,_,rt=arg['titles'][5]
            else:
                if multi[0]:
                    t,_,rt,_=arg['titles'][7]
                if multi[1]:
                    t,_,_,rt=arg['titles'][7]
        if g_type == 'SNPidx2':
            ylab = 'SNP-index2'
            key2 = 7
            if arg['--moving-avg'] != False:
                key = 11
            if arg['--distance-avg'] != False:
                key = 10
            if not multi[0] and not multi[1]:
                t,_,rt=arg['titles'][6]
            else:
                if multi[0]:
                    t,_,rt,_=arg['titles'][8]
                if multi[1]:
                    t,_,_,rt=arg['titles'][8]
    
    if g_type == 'DELTA':
        if not multi[0] and not multi[1]:
            t,rt=arg['titles'][9]
        else:
            if multi[0]:
                t,rt,_=arg['titles'][10]
            if multi[1]:
                t,_,rt=arg['titles'][10]
        if arg['--ref-genotype'] == True:
            ylab = '$\Delta$'+' (SNP-idx)'
            ticks_y=[-1, 0, 1]
            lim_y=(-1.2, 1.2)
            if arg['type'] == 'mbs':
                key2 = 16
            else:
                key2 = 18
        else:
            ylab = '|$\Delta$'+' (SNP-idx)|'
            ticks_y=[0, 0.25, 0.5, 0.75, 1]
            lim_y=(0, 1.2)
            if arg['type'] == 'mbs':
                key2 = 17
            else:
                key2 = 19
        if arg['--moving-avg'] != False:
            key = 20
        if arg['--distance-avg'] != False:
            key = 21
    return arg, ticks_y, lim_y, ylab, t, rt, key2, key

def AF_mono_graph(df, arg, g_type):
    sets = arg['sets']['mono']
    plt.rcParams['savefig.transparent'] = True
    chrom=arg['--chromosomes']
    typ=arg['--output-type']
    ticks_y=[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    ylab = 'SNP-index'
    boost = False
    multi = (False, False)
    arg, ticks_y, lim_y, ylab, t, rt, key2, key = chooseOptionsAFPlots(arg, multi, g_type, ticks_y, lim_y, ylab)

    for i in range(len(chrom)):
        
        d=df[df['#CHROM'] == chrom[i]]
        max_x=int(arg['contigs'][chrom[i]])
        x=d[['POS']]
        y=d[[g_type]]
        fig, ax=plt.subplots(figsize=(sets['figwidth'], sets['figheight']))
        ax.scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax.set(xlim=(0, max_x), ylim=lim_y)
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=sets['xticksSIZE'])
        ax.tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Genomic position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax.set_ylabel(ylabel=ylab, fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
        ax.set_yticks(ticks=ticks_y)
        ax.set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])
        ax.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax.spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax, g_type)
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax, g_type)
        if arg['type'] == 'mbs':
            if arg['--boost'] == True and arg['--moving-avg'] != False:
                keyB = 36
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plot_boost(d, arg, ax, max_x, sets, ticks_y, g_type, multi)
                    boost = True
            if arg['--boost'] == True and arg['--distance-avg'] != False:
                keyB = 35
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plotDistanceBoost(d, arg, ax, max_x, chrom[i], sets, ticks_y, g_type, multi)
                    boost = True
        
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
        plt.tight_layout()
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
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
            if arg['type'] == 'mbs':
                if boost == True:
                    if arg['--moving-avg'] != False:
                        cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--moving-avg'])))
                    if arg['--distance-avg'] != False:
                        cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--distance-avg'])))
            if arg['--ci95'] and g_type == 'DELTA' and (arg['--moving-avg'] != False or arg['--distance-avg'] != False):
                cap.append(arg['lines'][22].format(str(arg['n_markers'])))
            write_caption(f, cap, arg)
        del d

def pval_multi_Vertical_graph(df, arg):
    sets = arg['sets']['vertical']
    typ=arg['--output-type']
    min_y = min(df['log10PVALUE'])*1.05
    max_x = arg['max_lenght']
    t,rt,_ = arg['titles'][12]
    labs_list = list()
    f_name = list()
    threshold = log(0.05/len(df.axes[0]))/log(10)
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig,ax=plt.subplots(len(chrom), 1, figsize=(sets['figwidth'], sets['chromHeight']*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(sets['figwidth'], sets['chromHeight']*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['log10PVALUE']]
            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
            ax[i].set(xlim=(0, max_x), ylim=((min(min_y, threshold)-1)//1, 0))
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(labelbottom=False)
            ax[i].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            ax[i].set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
            ax[i].set_title('{})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'], length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            for spine in ax[i].spines.values():
                spine.set_linewidth(sets['axSpinesWIDTH'])
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], 'log10PVALUE')
            if arg['--distance-avg'] != False:
                plotDistanceAvg(d, chrom[i], arg, ax[i], 'log10PVALUE')
            if arg['--bonferroni'] and not d.empty:
                max_x_ch=(max(d['POS'])/max_x)
                ax[i].axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            if len(chrom) == 1:
                ax[1].remove()
            if i == len(chrom)-1 or len(chrom) == 1:
                ax[i].tick_params(labelbottom=True)
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=sets['xticksSIZE'], color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
                ax[i].set_xlabel(xlabel='Genomic position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])

        fig.tight_layout(pad=sets['pad'])
        fig.subplots_adjust(hspace=sets['hspace'])
        filename= rt + typ
        filename=check_save(arg, filename)
        f_name.append(filename)
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        cap = cap + labs_list
        cap.append(arg['lines'][12].format(arg['color_names']['dots'],''))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][14].format(arg['color_names']['log10PVALUE'],'',str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
                cap.append(arg['lines'][13].format(arg['color_names']['log10PVALUE'],'',str(arg['--distance-avg'])))
        if arg['--bonferroni'] and not d.empty:
            cap.append(arg['lines'][15].format(str(arg['n_markers'])))
        write_caption(f,cap,arg)
    del d


def ED_multi_Vertical_graph(df, arg):
    sets = arg['sets']['vertical']
    typ=arg['--output-type']
    ed100 = arg['ED100']
    if isinstance(ed100,pd.DataFrame):    
        max_y = max(ed100['ED100_4'])
    max_x = arg['max_lenght']
    t,rt,_ = arg['titles'][14]
    labs_list = list()
    f_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig,ax=plt.subplots(len(chrom), 1, figsize=(sets['figwidth'], sets['chromHeight']*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(sets['figwidth'], sets['chromHeight']*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['ED']]
            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
            ax[i].set(xlim=(0, max_x), ylim=(0, 1.5))
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='EDm', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
            ax[i].set_yticks(ticks=[0,0.5,1,1.5])
            ax[i].set_yticklabels(labels=[0,0.5,1,1.5], fontsize=sets['yticksSIZE'])
            ax[i].set_title('{})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
            ax[i].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            for spine in ax[i].spines.values():
                spine.set_linewidth(sets['axSpinesWIDTH'])
            ax[i].spines['top'].set_visible(False)
            if isinstance(ed100,pd.DataFrame):
                plot_ED100_4(arg, ax[i], max_x, max_y, ed100, chrom[i], sets)
            if len(chrom) == 1:
                ax[1].remove()
            if i == len(chrom)-1 or len(chrom) == 1:
                ax[i].tick_params(labelbottom=True)
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=sets['xticksSIZE'], color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
                ax[i].set_xlabel(xlabel='Genomic position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        fig.tight_layout(pad=sets['pad'])
        fig.subplots_adjust(hspace=sets['hspace'])
        filename= rt + typ
        filename=check_save(arg, filename)
        f_name.append(filename)
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
        del d
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        cap = cap + labs_list
        cap.append(arg['lines'][23].format(arg['color_names']['dots']))
        if isinstance(ed100,pd.DataFrame):
            cap.append(arg['lines'][24].format(arg['color_names']['mvg']))
        write_caption(f,cap,arg)

def EDmanhattanPlot(df, arg):
    sets = arg['sets']['manhattan']
    typ=arg['--output-type']
    ed100 = arg['ED100']
    if isinstance(ed100,pd.DataFrame):    
        max_y = max(ed100['ED100_4'])
    t,_,rt = arg['titles'][14]
    #labs_list = list()
    f_name = list()
    chrom = arg['--chromosomes']
    if len(arg['--chromosomes'])*sets['chromMinWIDTH'] <= sets['maxWIDTH']:
        af = len(arg['--chromosomes'])*sets['chromMinWIDTH']
    else:
        af = sets['maxWIDTH']
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, sets['figheight']))
    for i in range(len(chrom)):
        d=df[df['#CHROM'] == chrom[i]]
        max_x = int(arg['contigs'][chrom[i]])
        x=d[['POS']]
        y=d[['ED']]
        ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[i].set(xlim=(0, max_x), ylim=(0, 1.5))
        ax[i].set_xticks([])
        ax[i].set_xlabel(xlabel=f'{chrom[i]}', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax[i].set_ylabel(ylabel='EDm', fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'], rotation=90)
        #labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
        ax[i].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[i].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[i].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if isinstance(ed100,pd.DataFrame):
            flag = chrom[i] in ed100['#CHROM'].values
            ax2 = ax[i].twinx()
            c_ = arg['--palette']['mvg']
            ax2.set(xlim=(0, max_x), ylim=(0, 1.5))
            ax2.tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
            ax2.spines['top'].set_visible(False)
            ax2.axes.get_yaxis().set_visible(False)
            if flag:
                ed100ch = ed100[ed100['#CHROM'] == chrom[i]]
                ax2.plot(ed100ch['POS'], ed100ch['ED100_4'], c=c_, lw=arg['LINE_W'])
            if chrom[i] == chrom[-1]:
                ax2.set_yticks(ticks=[0, 0.5, 1, 1.5])
                ax2.set_yticklabels(labels=[0, 0.5, 1, 1.5], fontsize=sets['yticksSIZE'])
                ax2.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
                for spine in ax2.spines.values():
                    spine.set_linewidth(sets['axSpinesWIDTH'])
                ax2.axes.get_yaxis().set_visible(True)
                ax2.set_ylabel(ylabel='ED $\mathregular{100^4}$  $\mathregular{x10^8}$ ', fontsize=sets['xlabSIZE'], rotation=90, labelpad=sets['xlabSIZE'])
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        if len(chrom) == 1:
            ax[1].remove()
        if i > 0:
            ax[i].axes.get_yaxis().set_visible(False)
        ax[0].set_yticks(ticks=[0, 0.5, 1, 1.5])
        ax[0].set_yticklabels(labels=[0, 0.5, 1, 1.5], fontsize=sets['yticksSIZE'])
        del d
    #plt.tight_layout()
    fig.subplots_adjust(wspace=sets['wspace'], bottom=sets['bottom'])
    filename= rt + typ
    filename=check_save(arg, filename)
    f_name.append(filename)
    plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
    plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        #cap = cap + labs_list
        cap.append(arg['lines'][23].format(arg['color_names']['dots']))
        if isinstance(ed100,pd.DataFrame):
            cap.append(arg['lines'][24].format(arg['color_names']['mvg']))
        write_caption(f,cap,arg)

def G_multi_Vertical_graph(df, arg):
    sets = arg['sets']['vertical']
    typ=arg['--output-type']
    min_y, max_y = min(df['G']), max(df['G'])
    max_x =arg['max_lenght']
    t,rt,_ = arg['titles'][16]
    labs_list = list()
    f_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig,ax=plt.subplots(len(chrom), 1, figsize=(sets['figwidth'], sets['chromHeight']*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(sets['figwidth'], sets['chromHeight']*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['G']]
            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
            ax[i].set(xlim=(0, max_x), ylim=(min_y, max_y))
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(labelbottom=False)
            ax[i].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            ax[i].set_ylabel(ylabel='G-statistic', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
            ax[i].set_title('{})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'], length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            for spine in ax[i].spines.values():
                spine.set_linewidth(sets['axSpinesWIDTH'])
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
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=sets['xticksSIZE'], color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
                ax[i].set_xlabel(xlabel='Genomic position (bp)',  fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        fig.tight_layout(pad=sets['pad'])
        fig.subplots_adjust(hspace=sets['hspace'])
        filename= rt + typ
        filename=check_save(arg, filename)
        f_name.append(filename)
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
        del d
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        cap = cap + labs_list
        cap.append(arg['lines'][25].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][27].format(arg['color_names']['mvg'], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
                cap.append(arg['lines'][26].format(arg['color_names']['mvg'], str(arg['--distance-avg'])))
        write_caption(f,cap,arg)

def GmanhattanPlot(df, arg):
    sets = arg['sets']['manhattan']
    typ=arg['--output-type']
    min_y, max_y = min(df['G']), max(df['G'])
    t,_,rt = arg['titles'][16]
    #labs_list = list()
    f_name = list()
    chrom = arg['--chromosomes']
    if len(arg['--chromosomes'])*sets['chromMinWIDTH'] <= sets['maxWIDTH']:
        af = len(arg['--chromosomes'])*sets['chromMinWIDTH']
    else:
        af = sets['maxWIDTH']
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, sets['figheight']))
    for i in range(len(chrom)):
        d=df[df['#CHROM'] == chrom[i]]
        max_x = int(arg['contigs'][chrom[i]])
        x=d[['POS']]
        y=d[['G']]
        ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[i].set(xlim=(0, max_x), ylim=(min_y, max_y))
        ax[i].set_xticks([])
        ax[i].set_xlabel(xlabel=f'{chrom[i]}', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax[i].set_ylabel(ylabel='G-statistic', fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'], rotation=90)
        #labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
        ax[i].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[i].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[i].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[i], 'G')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[i], 'G')

        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        if len(chrom) == 1:
            ax[1].remove()
        if i > 0:
            ax[i].axes.get_yaxis().set_visible(False)
        del d
    #fig.tight_layout()
    fig.subplots_adjust(wspace=sets['wspace'], bottom=sets['bottom'])
    filename= rt + typ
    filename=check_save(arg, filename)
    f_name.append(filename)
    plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
    plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        #cap = cap + labs_list
        cap.append(arg['lines'][25].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][27].format(arg['color_names']['mvg'], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][26].format(arg['color_names']['mvg'], str(arg['--distance-avg'])))
        write_caption(f,cap,arg)

def AF_multi_Vertical_graph(df, arg, g_type):
    sets = arg['sets']['vertical']
    typ=arg['--output-type']
    max_x=arg['max_lenght']
    ticks_y=[0, 0.5, 1]#[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    ylab = 'SNP-index'
    boost = False
    multi = (True,False)    #Multivertical Not Manhattan
    labs_list = list()
    f_name = list()
    arg, ticks_y, lim_y, ylab, t, rt, key2, key = chooseOptionsAFPlots(arg, multi, g_type, ticks_y, lim_y, ylab)

    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig, ax=plt.subplots(len(chrom), 1, figsize=(sets['figwidth'], sets['chromHeight']*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(sets['figwidth'],sets['chromHeight']*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[[g_type]]

            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
            ax[i].set(xlim=(0, max_x), ylim=lim_y)
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel=ylab, fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
            ax[i].set_title('{})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x = sets['titleXPOS'], y=sets['titleYPOS'])
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].set_yticks(ticks=ticks_y)
            ax[i].set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])
            ax[i].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            for spine in ax[i].spines.values():
                spine.set_linewidth(sets['axSpinesWIDTH'])
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], g_type)
            if arg['--distance-avg'] != False:
                plotDistanceAvg(d, chrom[i], arg, ax[i], g_type)
            if arg['type'] == 'mbs':
                if arg['--boost'] == True and arg['--moving-avg'] != False:
                    keyB = 36
                    if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                        plot_boost(d, arg, ax[i], max_x, sets, ticks_y, g_type,multi)
                        boost = True
                if arg['--boost'] == True and arg['--distance-avg'] != False:
                    keyB = 35
                    if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                        plotDistanceBoost(d, arg, ax[i], max_x, chrom[i], sets, ticks_y, g_type,multi)
                        boost = True
            
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
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=sets['xticksSIZE'])
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
                ax[i].set_xlabel(xlabel='Genomic position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
            del d
        fig.tight_layout(pad=sets['pad'])
        fig.subplots_adjust(hspace=sets['hspace'])
        filename= rt + typ
        filename=check_save(arg, filename)
        f_name.append(filename)
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
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
        if arg['type'] == 'mbs':
            if boost == True:
                if arg['--moving-avg'] != False:
                    cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--moving-avg'])))
                if arg['--distance-avg'] != False:
                    cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--distance-avg'])))
        if arg['--ci95'] and g_type == 'DELTA' and (arg['--moving-avg'] != False or arg['--distance-avg'] != False):
            cap.append(arg['lines'][22].format(str(arg['n_markers'])))
        write_caption(f, cap,arg)

def AF12_multi_Vertical_graph(df, arg):
    sets = arg['sets']['vertical']
    typ=arg['--output-type']
    max_x=arg['max_lenght']
    if arg['type'] == 'mbs':
        k1, k2 = 0,1
        ti1, ti2 = 2,3
    else:
        k1,k2 = 6,7
        ti1,ti2 = 7,8
    ticks_y=[0, 0.5, 1]#[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    t1,rt1,_,_=arg['titles'][ti1]
    t2,rt2,_,_=arg['titles'][ti2]
    labs_list = list()
    f1_name = list()
    f2_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig, ax=plt.subplots(len(chrom), 1, figsize=(sets['figwidth'], sets['chromHeight']*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(sets['figwidth'], sets['chromHeight']*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['SNPidx1']]
            #AF1&2
            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
            ax[i].set(xlim=(0, max_x), ylim=lim_y)
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='SNP-index1', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x =sets['titleXPOS'], y=sets['titleYPOS']) 
            labs_list.append('{}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))  
            ax[i].set_yticks(ticks=ticks_y)
            ax[i].set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])
            ax[i].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            for spine in ax[i].spines.values():
                spine.set_linewidth(sets['axSpinesWIDTH'])
            if arg['--moving-avg'] != False:
                if arg['type'] == 'mbs':
                    l1,l2 = 3,5
                else:
                    l1,l2 = 9,11
                plot_avg(d, arg, ax[i], 'SNPidx2')
                plot_avg(d, arg, ax[i], 'SNPidx1')
            if arg['--distance-avg'] != False:
                if arg['type'] == 'mbs':
                    l1,l2 = 2,4
                else:
                    l1,l2 = 8,10
                plotDistanceAvg(d, chrom[i], arg, ax[i], 'SNPidx2')
                plotDistanceAvg(d, chrom[i], arg, ax[i], 'SNPidx1')
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            if len(chrom) == 1:
                ax[1].remove()
            if i == len(chrom)-1 or len(chrom) == 1:
                ax[i].tick_params(labelbottom=True)
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=sets['xticksSIZE'], color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
                ax[i].set_xlabel(xlabel='Genomic position (bp)', fontsize=sets['xlabSIZE'],labelpad=sets['xlabDIST'])
            del d
        plt.tight_layout(pad=sets['pad'])
        fig.subplots_adjust(hspace=sets['hspace'])
        filename=rt1+typ
        filename=check_save(arg, filename)
        f1_name.append(filename)
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt1)
        cap.append(', '.join(f1_name))
        cap.append(t1)
        cap = cap + labs_list
        cap.append(arg['lines'][k1].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][l1].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])))
            cap.append(arg['lines'][l2].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][l1].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])))
            cap.append(arg['lines'][l2].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
        write_caption(f,cap,arg)
    
    #AF1&2
    f2_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig,ax=plt.subplots(len(chrom), 1, figsize=(sets['figwidth'], sets['chromHeight']*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(sets['figwidth'], sets['chromHeight']*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['SNPidx2']]
            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
            ax[i].set(xlim=(0, max_x), ylim=lim_y)
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='SNP-index2', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x =sets['titleXPOS'], y=sets['titleYPOS'])
            ax[i].set_yticks(ticks=ticks_y)
            ax[i].set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])
            ax[i].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            for spine in ax[i].spines.values():
                spine.set_linewidth(sets['axSpinesWIDTH'])
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
                ax[i].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=sets['xticksSIZE'], color='black')
                ax[i].xaxis.set_major_formatter(ScalarFormatter())
                ax[i].ticklabel_format(axis='x', style='scientific', scilimits=(6, 6), useMathText=True)
                ax[i].xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
                ax[i].set_xlabel(xlabel='Genomic position (bp)', fontsize=sets['xlabSIZE'],labelpad=sets['xlabDIST'])

        fig.tight_layout(pad=sets['pad'])
        fig.subplots_adjust(hspace=sets['hspace'])
        filename=rt2+typ
        filename=check_save(arg, filename)
        f2_name.append(filename)
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()

    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt2)
        cap.append(', '.join(f2_name))
        cap.append(t2)
        cap = cap + labs_list
        cap.append(arg['lines'][k2].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][l2].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])))
            cap.append(arg['lines'][l1].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
            cap.append(arg['lines'][l2].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])))
            cap.append(arg['lines'][l1].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])))
        write_caption(f,cap,arg)


def combinedPlot(df, arg):
    sets = arg['sets']['combined']
    global fields
    chrom=arg['--chromosomes']
    typ=arg['--output-type']
    t,rt = arg['titles'][17]
    threshold=log(0.05/len(df.axes[0]))/log(10)
    min_yp=min(df['log10PVALUE'])*1.05
    if arg['type'] == 'mbs':
        key = 30
    else:
        key = 31
    if arg['--moving-avg'] != False:
            k1 = 29
    if arg['--distance-avg'] != False:
            k1 = 28
    for i in range(len(chrom)):
        fig, ax=plt.subplots(3, 2, figsize=(sets['figwidth'], sets['figheight']))
        d=df[df['#CHROM'] == chrom[i]]
        max_x = arg['contigs'][chrom[i]]
        if max_x > 60000000:
            paso = 10e6
        else:
            paso = 5e6
        x = d[['POS']]
        y1 = d['SNPidx1']
        y2 = d['SNPidx2']
        y3 = d['DELTA']
        y4 = d['ED']
        y5 = d['G']
        y6 = d['log10PVALUE']

        #SNPidx1
        ax[0,0].scatter(x, y1, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[0,0].set(xlim=(0, max_x), ylim=(0, 1))
        ax[0,0].set_xticks(ticks=np.arange(0, max_x, paso))
        ax[0,0].tick_params(labelbottom=False)
        ax[0,0].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax[0,0].set_yticks(ticks=[0, 0.5, 1])
        ax[0,0].set_yticklabels(labels=[0, 0.5, 1], fontsize=sets['yticksSIZE'])
        ax[0,0].set_ylabel('SNP-index1', fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'])
        ax[0,0].yaxis.set_label_coords(sets['ylabXPOS'],sets['ylabYPOS'])
        ax[0,0].set_title('a)', fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
        ax[0,0].axhline(y=0.5, color='black', linestyle='dashed', linewidth=0.75)
        ax[0,0].spines['top'].set_visible(False)
        ax[0,0].spines['right'].set_visible(False)
        ax[0,0].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[0,0].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[0,0].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[0,0], 'SNPidx1')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[0,0], 'SNPidx1')

        # SNPidx2
        ax[1,0].scatter(x, y2, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[1,0].set(xlim=(0, max_x), ylim=(0, 1))
        ax[1,0].set_xticks(ticks=np.arange(0, max_x, paso))
        ax[1,0].tick_params(labelbottom=False)
        ax[1,0].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax[1,0].set_yticks(ticks=[0, 0.5, 1])
        ax[1,0].set_yticklabels(labels=[0, 0.5, 1], fontsize=sets['yticksSIZE'])
        ax[1,0].set_ylabel('SNP-index2', fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'])
        ax[1,0].yaxis.set_label_coords(sets['ylabXPOS'],sets['ylabYPOS'])
        ax[1,0].set_title('b)', fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
        ax[1,0].axhline(y=0.5, color='black', linestyle='dashed', linewidth=0.75)
        ax[1,0].spines['top'].set_visible(False)
        ax[1,0].spines['right'].set_visible(False)
        ax[1,0].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[1,0].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[1,0].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[1,0], 'SNPidx2')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[1,0], 'SNPidx2')

        # DELTA
        if arg['--ref-genotype'] == True:
            ylab = '$\Delta$'+' (SNP-index)'
            ticks_y=[-1, 0, 1]
            lim_y=(-1.2, 1.2)
        else:
            ylab = '|$\Delta$'+' (SNP-index)|'
            ticks_y=[0, 0.25, 0.5, 0.75, 1]
            lim_y=(0, 1.2)
        ax[2,0].scatter(x, y3, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[2,0].set(xlim=(0, max_x), ylim=lim_y)
        ax[2,0].set_xticks(ticks=np.arange(0, max_x, paso))
        ax[2,0].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax[2,0].set_xticklabels(labels=np.arange(0, max_x, paso),fontsize=sets['xticksSIZE'])
        ax[2,0].set_ylabel(ylabel= ylab ,fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'])
        ax[2,0].yaxis.set_label_coords(sets['ylabXPOS'],sets['ylabYPOS'])
        ax[2,0].set_title('c)', fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
        ax[2,0].set_yticks(ticks=ticks_y)
        ax[2,0].set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])

        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[2,0], 'DELTA')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[2,0], 'DELTA')
        if arg['--ci95'] and (arg['--moving-avg'] != False or arg['--distance-avg'] != False):
            if arg['--moving-avg'] != False:
                calc_ci(d, arg, ax[2,0])
            if arg['--distance-avg'] != False:
                distanceCI(d, arg, ax[2,0], chrom[i])
        ax[2,0].axhline(y=0, color = 'black', linestyle='dashed', linewidth=0.75)
        ax[2,0].xaxis.set_major_formatter(ScalarFormatter())
        ax[2,0].ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax[2,0].xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax[2,0].set_xlabel(xlabel='Genomic position (bp)',fontsize=sets['xlabSIZE'],labelpad=sets['xlabDIST'])
        ax[2,0].spines['top'].set_visible(False)
        ax[2,0].spines['right'].set_visible(False)
        ax[2,0].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[2,0].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[2,0].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])

        # ED
        ed100 = arg['ED100']
        if isinstance(ed100,pd.DataFrame):
            max_y = max(ed100['ED100_4'])
        ax[0,1].scatter(x, y4, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[0,1].set(xlim=(0, max_x), ylim=(0, 1.5))
        ax[0,1].set_xticks(ticks=np.arange(0, max_x, paso))
        ax[0,1].tick_params(labelbottom=False)
        ax[0,1].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax[0,1].set_ylabel(ylabel='EDm', fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'], rotation=90)
        ax[0,1].yaxis.set_label_coords(sets['ylabXPOS'],sets['ylabYPOS'])
        ax[0,1].set_title('d)',fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
        ax[0,1].set_yticks(ticks=[0, 0.5, 1, 1.5])
        ax[0,1].set_yticklabels(labels=[0, 0.5, 1, 1.5], fontsize=sets['yticksSIZE'])
        ax[0,1].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[0,1].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[0,1].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if isinstance(ed100,pd.DataFrame):
            plot_ED100_4(arg, ax[0,1], max_x, max_y, ed100, chrom[i], sets)
        ax[0,1].spines['top'].set_visible(False)
        ax[0,1].spines['right'].set_visible(False)

        # G-statistic
        ax[1,1].scatter(x, y5, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[1,1].set(xlim=(0, max_x), ylim=(min(df['G']), max(df['G'])))
        ax[1,1].set_xticks(ticks=np.arange(0, max_x, paso))
        ax[1,1].tick_params(labelbottom=False)
        ax[1,1].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax[1,1].set_ylabel(ylabel='G-statistic', fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'], rotation=90)
        ax[1,1].yaxis.set_label_coords(sets['ylabXPOS'],sets['ylabYPOS'])
        ax[1,1].set_title('e)', fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
        ax[1,1].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[1,1].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[1,1].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[1,1], 'G' )
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[1,1], 'G')
        ax[1,1].spines['top'].set_visible(False)
        ax[1,1].spines['right'].set_visible(False)


        # log10PVALUE
        ax[2,1].scatter(x, y6, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[2,1].set(xlim=(0, max_x), ylim=((min(min_yp, threshold)-1)//1, 0))
        ax[2,1].set_xticks(ticks=np.arange(0, max_x, paso))
        ax[2,1].set_xticklabels(labels=np.arange(0, max_x, paso),fontsize=sets['xticksSIZE'])
        ax[2,1].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax[2,1].set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)',fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'])
        ax[2,1].yaxis.set_label_coords(sets['ylabXPOS'],sets['ylabYPOS'])
        ax[2,1].set_title('f)', fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
        ax[2,1].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[2,1], 'log10PVALUE')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[2,1], 'log10PVALUE')
        if arg['--bonferroni'] and not d.empty:
            max_x_ch=(max(d['POS'])/max_x)
            ax[2,1].axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        ax[2,1].xaxis.set_major_formatter(ScalarFormatter())
        ax[2,1].ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax[2,1].xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax[2,1].set_xlabel(xlabel='Genomic position (bp)',fontsize=sets['xlabSIZE'],labelpad=sets['xlabDIST'])
        ax[2,1].spines['top'].set_visible(False)
        ax[2,1].spines['right'].set_visible(False)
        ax[2,1].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[2,1].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[2,1].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        # Save file
        fig.tight_layout()
        fig.subplots_adjust(hspace=sets['hspace'], wspace=sets['wspace'])
        rtch = rt.format(chrom[i])
        filename=rtch + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
        cap = list()

        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t.format(chrom[i]))
            if arg['--moving-avg'] != False or arg['--distance-avg'] != False:  
                cap.append(arg['lines'][key].format(
                                        arg['lines'][k1].format(str(arg['--moving-avg'])),\
                                        '' if arg['--ref-genotype'] == True else ', in absolute value',\
                                        arg['lines'][22].format(str(arg['n_markers'])) if arg['--ci95'] else '',\
                                        arg['lines'][24].format(arg['color_names']['mvg']) if isinstance(ed100,pd.DataFrame) else '',\
                                        arg['lines'][15].format(str(arg['n_markers'])) if arg['--bonferroni'] else ''
                                        ))
            else:
                cap.append(arg['lines'][key].format(
                                        '',\
                                        '' if arg['--ref-genotype'] == True else ', in absolute value',\
                                        '',\
                                        arg['lines'][24].format(arg['color_names']['mvg']) if isinstance(ed100,pd.DataFrame) else '',\
                                        arg['lines'][15].format(str(arg['n_markers'])) if arg['--bonferroni'] else ''
                                        ))
            write_caption(f,cap, arg)
            del d


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
    
    ax.plot(d['avgx'], d['avgy'], c=c_, lw=arg['LINE_W'])#lw=0.75

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

    ax.plot(res['POSx'], res['VALy'], c=c_, lw=arg['LINE_W'])


def plot_boostManhattan(d, arg, ax, max_x, chrom, i, sets, ticks_y, g_type):
    # Mediamovil Boost
    d['DP']=d['DPdom_2']+d['DPrec_2']
    d['prody']=d['SNPidx2'] * d['DP']
    d['medboost']=(d['prody'].rolling(arg['--moving-avg']).sum() / \
                   d['DP'].rolling(arg['--moving-avg']).sum())
    d['BOOST2'] = 1/( arg['lim']+abs(1-1/(np.maximum(d['medboost'], 1-d['medboost']))))
    d['BOOST2'] = d['BOOST2']* arg['lim']
    d['prodboostx']=d['POS'] * d['DP']
    d['medboostx']=(d['prodboostx'].rolling(
        arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())
    c_=arg['--palette']['BOOST']
    if g_type == 'MAX_SNPidx2':
        ticks_y = [0,0.25,0.5,0.75,1]
    ax2=ax.twinx()
    ax2.spines['top'].set_visible(False)
    ax2.set(xlim=(0, max_x), ylim=(0, 1))
    ax2.set_yticks([])
    ax2.plot(d['medboostx'], d['BOOST2'], c=c_, lw=arg['LINE_W'], linestyle='dashed')
    ax2.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
    for spine in ax2.spines.values():
        spine.set_linewidth(sets['axSpinesWIDTH'])

    if chrom[i] == chrom[-1]:
        ax2.axes.get_yaxis().set_visible(True)
        ax2.set_yticks(ticks=ticks_y)
        ax2.set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])
        ax2.set_ylabel(ylabel='Boost',fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])

def plot_boost(d, arg, ax, max_x, sets, ticks_y, g_type, multi):
    # Mediamovil Boost
    d['DP']=d['DPdom_2']+d['DPrec_2']
    d['prody']=d['SNPidx2'] * d['DP']
    d['medboost']=(d['prody'].rolling(arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())
    d['BOOST2'] = 1/( arg['lim']+abs(1-1/(np.maximum(d['medboost'], 1-d['medboost']))))
    d['BOOST2'] = d['BOOST2']* arg['lim']
    d['prodboostx']=d['POS'] * d['DP']
    d['medboostx']=(d['prodboostx'].rolling(
        arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())
    ax2=ax.twinx()
    if g_type == 'MAX_SNPidx2':
        ticks_y =  [0,0.5,1] if multi[0] else [0,0.25,0.5,0.75,1]
    ax2.spines['top'].set_visible(False)
    c_=arg['--palette']['BOOST']
    ax2.plot(d['medboostx'], d['BOOST2'], c=c_, lw=arg['LINE_W'], linestyle='dashed')
    ax2.set(xlim=(0, max_x), ylim=(0, 1))
    ax2.set_yticks(ticks=ticks_y)
    ax2.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
    for spine in ax2.spines.values():
        spine.set_linewidth(sets['axSpinesWIDTH'])
    ax2.set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])
    ax2.set_ylabel(ylabel='Boost', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])

def plotDistanceBoostManhattan(d, arg, ax, max_x, chrom, i, sets, ticks_y, g_type):
    # Mediamovil Boost
    d['DP']=d['DPdom_2']+d['DPrec_2']
    d['prodx'] = d['POS']*d['DP']
    d['prody'] = d['SNPidx2']*d['DP']
    interval_ranges = [x for x in range(0, max_x, arg['--distance-avg'])]
    if interval_ranges[-1] < max_x:
        interval_ranges.append(max_x)
    interval_labels = [f'{start}-{end-1}' for start, end in zip(interval_ranges, interval_ranges[1:])]
    d['Interval'] = pd.cut(d['POS'], bins=interval_ranges, labels=interval_labels)
    res = d.groupby('Interval').agg({'prodx': sum, 'prody': sum, 'DP': sum})
    res['BOOST2'] = 1/( arg['lim'] + abs(1-1/(np.maximum(res['prody'], 1-res['prody']))))
    res['BOOST2'] = res['BOOST2']*arg['lim']
    res['POSx'] = (res['prodx']/res['DP'])
    res['VALy'] = (res['BOOST2']/res['DP'])
    res = res.dropna()

    c_=arg['--palette']['BOOST']
    ax2=ax.twinx()
    ax2.spines['top'].set_visible(False)
    ax2.set(xlim=(0, max_x), ylim=(0, 1))
    ax2.set_yticks([])
    ax2.plot(res['POSx'], res['VALy'], c=c_, lw=arg['LINE_W'], linestyle='dashed')
    ax2.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
    for spine in ax2.spines.values():
        spine.set_linewidth(sets['axSpinesWIDTH'])
    if g_type == 'MAX_SNPidx2':
        ticks_y = [0,0.25,0.5,0.75,1]
    if chrom[i] == chrom[-1]:
        ax2.axes.get_yaxis().set_visible(True)
        ax2.set_yticks(ticks=ticks_y)
        ax2.set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])
        ax2.set_ylabel(ylabel='Boost',fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])


def plotDistanceBoost(d, arg, ax, max_x, chrom, sets, ticks_y, g_type, multi):
    # Mediamovil Boost
    max_x = arg['contigs'][chrom]
    d['DP']=d['DPdom_2']+d['DPrec_2']
    d['prodx'] = d['DP']*d['POS']
    d['prody'] = d['DP']*d['SNPidx2']
    interval_ranges = [x for x in range(0, max_x, arg['--distance-avg'])]
    if interval_ranges[-1] < max_x:
        interval_ranges.append(max_x)
    interval_labels = [f'{start}-{end-1}' for start, end in zip(interval_ranges, interval_ranges[1:])]
    d['Interval'] = pd.cut(d['POS'], bins=interval_ranges, labels=interval_labels)
    res = d.groupby('Interval').agg({'prodx': sum, 'prody': sum, 'DP': sum})
    res['BOOST2'] = 1/( arg['lim'] + abs(1-1/(np.maximum(res['prody'], 1-res['prody']))))
    res['BOOST2'] = res['BOOST2']*arg['lim']
    res['POSx'] = (res['prodx']/res['DP'])
    res['VALy'] = (res['BOOST2']/res['DP'])
    res = res.dropna()
    ax2=ax.twinx()
    if g_type == 'MAX_SNPidx2':
        ticks_y =  [0,0.5,1] if multi[0] else [0,0.25,0.5,0.75,1]
    ax2.spines['top'].set_visible(False)
    c_=arg['--palette']['BOOST']
    ax2.plot(res['POSx'], res['VALy'], c=c_, lw=arg['LINE_W'], linestyle='dashed')
    ax2.set(xlim=(0, max_x), ylim=(0, 1))
    ax2.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
    for spine in ax2.spines.values():
        spine.set_linewidth(sets['axSpinesWIDTH'])
    ax2.set_yticks(ticks=ticks_y)
    ax2.set_yticklabels(labels=ticks_y, fontsize=sets['yticksSIZE'])
    ax2.set_ylabel(ylabel='Boost', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])

def plot_ED100_4(arg, ax, max_x, max_y, ed100, ch, sets):
    flag = ch in ed100['#CHROM'].values
    ax2 = ax.twinx()
    c_ = arg['--palette']['mvg']
    ax2.set(xlim=(0, max_x), ylim=(0, max_y))
    ax2.set_yticks(ticks=[0, 0.5, 1, 1.5])
    ax2.set_yticklabels(labels=[0, 0.5, 1, 1.5], fontsize=sets['yticksSIZE'])
    ax2.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
    for spine in ax2.spines.values():
        spine.set_linewidth(sets['axSpinesWIDTH'])
    if flag:
        ed100ch = ed100[ed100['#CHROM'] == ch]
        #print(ed100ch)
        ax2.set_ylabel(ylabel='ED$\mathregular{100^4}$$\mathregular{x10^8}$ ', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
        ax2.plot(ed100ch['POS'], ed100ch['ED100_4'], c=c_, lw=arg['LINE_W'])
    ax2.spines['top'].set_visible(False)

def create_caption(arg, res_tit):
    f= io.open(arg['captions_dir'] + res_tit + '.txt', 'a', encoding='utf8')
    return f

def create_subfolder(arg, sf):
    arg['sub_folder'] = os.path.abspath(arg['--outdir'])+'/'+sf
    try:
        os.makedirs(arg['sub_folder'])
    except FileExistsError:
        pass
    
    if arg['--captions']:
        arg['captions_dir'] = arg['sub_folder']+'captions/'
        try:
            os.makedirs(arg['captions_dir'])
        except FileExistsError:
            pass

    return arg

def write_caption(f, text, arg):
    for sentence in text:
        f.write(sentence+' ')
    if '--window' in arg.keys():
        f.write('The dots correspond to genomic regions defined by non-overlapping bins of {} consecutive markers. '.format(str(arg['--window'])))
    if '--distance' in arg.keys():
        f.write('The dots correspond to bins of all markers within a range of {} pb. '.format(str(arg['--distance'])))
    f.write('\n'*2)
    f.close()
