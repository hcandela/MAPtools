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
    read_palette(arg)
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
        #if 'SNPidx1' in arg['--fields'] and 'SNPidx2' in arg['--fields']:
        #    arg['--combine'] = True
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
                    print(f'Warning: there is not enough markers in {ch} to calculate ED100', file=sys.stderr)
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

def read_palette(arg):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    json_file_path = os.path.join(script_dir,'..','palette.json')
    try:
        with open(json_file_path) as json_file:
            data = json.load(json_file)
    except FileNotFoundError:
        print('Error: Please put the file \"palette.json\" in the MAPtools folder.', file=sys.stderr)
        sys.exit()
    #arg['DPI'] = data['DPI']
    #arg['DOT_SIZE'] = data['DOT_SIZE']
    #arg['LINE_W'] = data['LINE_W']
    palette = arg['--palette']
    if palette not in data.keys():
        print(f'Warning: {palette} is not one of the available palette names.\nSwitching to standard', file=sys.stderr)
        palette = 'standard'
    
    for key,val in data[palette].items():
        if val == list():
            palette = 'standard'
            break
    arg['--palette'] = {k: v[0] for k, v in data[palette].items()}
    arg['color_names'] = {k: v[1] for k, v in data[palette].items()}
    
def read_settings(arg):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    json_file_path = os.path.join(script_dir,'..','settings.json')
    try:
        with open(json_file_path) as json_file:
            data = json.load(json_file)
    except FileNotFoundError:
        print('Error: Please put the file \"settings.json\" in the MAPtools folder.', file=sys.stderr)
        sys.exit()

    for key,dicc in data.items():
        for k,v in dicc.items():
            if v == 'None':
                dicc[k] = None
    arg['sets'] = data

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
                if 'BOOST' in arg['--fields']:
                    boost = 1/(sys.float_info.min + abs(1 - 1/max(rSNPidx2, 1-rSNPidx2)))
                    res += [boost]

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
        fig, ax = plt.subplots(figsize=(10, 4.2))
        plt.subplots_adjust(bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
        ax.scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax.set(xlim=(0, max_x), ylim=(0, 1.5))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=sets['xticksSIZE'])
        ax.tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Chromosomal position (pb)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax.set_ylabel(ylabel='Euclidean distance', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
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
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()

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
        fig, ax = plt.subplots(figsize=(10, 4.2))
        plt.subplots_adjust(bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
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
            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
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


def AF_manhattan_plot(df, arg, g_type):
    sets = arg['sets']['manhattan']
    typ = arg['--output-type']
    #labs_list = list()
    f_name = list()
    ylab = 'Allele Frequency'
    boost = False
    ticks_y=[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    multi = (False, True)
    arg, ticks_y, lim_y, ylab, t, rt, key2, key = chooseOptionsAFPlots(arg, multi, g_type, ticks_y, lim_y, ylab)
    chrom = arg['--chromosomes']
    if len(arg['--chromosomes'])*sets['chromMinWIDTH'] <= sets['maxWIDTH']:
        af = len(arg['--chromosomes'])*sets['chromMinWIDTH']
    else:
        af = sets['maxWIDTH']
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, 3))

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
            if arg['--boost'] != False:
                keyB = 36
                if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                    plot_boostManhattan(d, arg, ax[i], max_x, chrom, i)
                    boost = True
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plot_boostManhattan(d, arg, ax[i], max_x, chrom, i)
                    boost = True
            if arg['--distance-boost'] != False:
                keyB = 35
                if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                    plotDistanceBoostManhattan(d, arg, ax[i], max_x, chrom, i)
                    boost = True
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plotDistanceBoostManhattan(d, arg, ax[i], max_x, chrom, i)
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
    fig.subplots_adjust(wspace=0, bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
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
                if arg['--boost'] != False:
                    cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--boost'])))
                if arg['--distance-boost'] != False:
                    cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--distance-boost'])))
        if arg['--ci95'] and g_type == 'DELTA':
                cap.append(arg['lines'][22].format(str(arg['n_markers'])))
        write_caption(f, cap,arg)


def AFCombinedManhattanPlot(df, arg):
    sets = arg['sets']['manhattan']
    typ = arg['--output-type']
    #labs_list = list()
    f_name = list()
    ylab = 'Allele Frequency'
    ticks_y=[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    t,rt = arg['titles'][4]

    chrom = arg['--chromosomes']
    if len(arg['--chromosomes'])*sets['chromMinWIDTH'] <= sets['maxWIDTH']:
        af = len(arg['--chromosomes'])*sets['chromMinWIDTH']
    else:
        af = sets['maxWIDTH']
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, 3))

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
    fig.subplots_adjust(wspace=0, bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
    filename = rt + typ
    filename = check_save(arg, filename)
    f_name.append(filename)
    plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
    plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        #cap.append(t.format(', '.join(labs_list)))
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
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, 3))
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
        if arg['--bonferroni']:
            max_x_ch = (max(d['POS'])/max_x)
            ax[i].axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        if len(chrom) == 1:
            ax[1].remove()
        if i > 0:
            ax[i].axes.get_yaxis().set_visible(False)
    fig.subplots_adjust(wspace=0, bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
    filename = rt + typ
    filename = check_save(arg, filename)
    print(filename)
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
        fig, ax=plt.subplots(figsize=(10, 4.2))
        plt.subplots_adjust(bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
        ax.scatter(x, y, s=arg['DOT_SIZE'], color=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax.set(xlim=(0, max_x), ylim=((min(min_y, threshold)-1)//1, 0))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=sets['xticksSIZE'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Chromosomal position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax.set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)', fontsize=sets['ylabSIZE'], labelpad=sets['ylabDIST'])
        ax.tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax.spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax, 'log10PVALUE')
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax, 'log10PVALUE')
        if arg['--bonferroni']:
            max_x_ch=(max(d['POS'])/max_x)
            ax.axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        rtch = rt.format(chrom[i])

        filename=rtch + typ
        filename=check_save(arg, filename)

        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()

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
        fig, ax=plt.subplots(figsize=(10, 4.2))
        plt.subplots_adjust(bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
        ax.scatter(x, y1, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax.set(xlim=(0, max_x), ylim=(0, 1))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=sets['xticksSIZE'])
        ax.tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Chromosomal position (bp)'.format(chrom[i]), fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax.set_ylabel(ylabel='Allele Frequency', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
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
        fig, ax=plt.subplots(figsize=(10, 4.2))
        ax.scatter(x, y2, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        plt.subplots_adjust(bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
        ax.set(xlim=(0, max_x), ylim=(0, 1))
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=sets['xticksSIZE'])
        ax.tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Chromosomal position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax.set_ylabel(ylabel='Allele frequency', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
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


def chooseOptionsAFPlots(arg, multi, g_type, ticks_y, lim_y, ylab):
    key = False
    if arg['type'] == 'mbs':
        if g_type == 'SNPidx1':
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
            ylab = 'Max A. F.'
            if not multi[0] and not multi[1]:
                t,rt=arg['titles'][18]
            else:
                if multi[0]:
                    t,rt,_=arg['titles'][19]
                if multi[1]:
                    t,_,rt=arg['titles'][19]
    else:
        if g_type == 'SNPidx1':
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
    ylab = 'Allele Frequency'
    boost = False
    multi = (False, False)
    arg, ticks_y, lim_y, ylab, t, rt, key2, key = chooseOptionsAFPlots(arg, multi, g_type, ticks_y, lim_y, ylab)

    for i in range(len(chrom)):
        
        d=df[df['#CHROM'] == chrom[i]]
        max_x=int(arg['contigs'][chrom[i]])
        x=d[['POS']]
        y=d[[g_type]]
        fig, ax=plt.subplots(figsize=(10, 4.2))
        plt.subplots_adjust(bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
        ax.scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax.set(xlim=(0, max_x), ylim=lim_y)
        ax.set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax.set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=sets['xticksSIZE'])
        ax.tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax.xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax.set_xlabel(xlabel='Chromosomal position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
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
            if arg['--boost'] != False:
                keyB = 36
                if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                    plot_boost(d, arg, ax, max_x, sets)
                    boost = True
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plot_boost(d, arg, ax, max_x, sets)
                    boost = True
            if arg['--distance-boost'] != False:
                keyB = 35
                if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                    plotDistanceBoost(d, arg, ax, max_x, chrom[i], sets)
                    boost = True
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plotDistanceBoost(d, arg, ax, max_x, chrom[i], sets)
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
                    if arg['--boost'] != False:
                        cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--boost'])))
                    if arg['--distance-boost'] != False:
                        cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--distance-boost'])))
            if arg['--ci95'] and g_type == 'DELTA':
                cap.append(arg['lines'][22].format(str(arg['n_markers'])))
            write_caption(f, cap, arg)

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
            fig,ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
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
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'], length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            for spine in ax[i].spines.values():
                spine.set_linewidth(sets['axSpinesWIDTH'])
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], 'log10PVALUE')
            if arg['--distance-avg'] != False:
                plotDistanceAvg(d, chrom[i], arg, ax[i], 'log10PVALUE')
            if arg['--bonferroni']:
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
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])


        fig.subplots_adjust(hspace=sets['hspace'], bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
        filename= rt + typ
        filename=check_save(arg, filename)
        f_name.append(filename)
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t.format(' '.join(labs_list)))
        cap = cap + labs_list
        cap.append(arg['lines'][12].format(arg['color_names']['dots'],''))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][14].format(arg['color_names']['log10PVALUE'],'',str(arg['--moving-avg'])))
        if arg['--distance-avg'] != False:
                cap.append(arg['lines'][13].format(arg['color_names']['log10PVALUE'],'',str(arg['--distance-avg'])))
        if arg['--bonferroni']:
            cap.append(arg['lines'][15].format(str(arg['n_markers'])))
        write_caption(f,cap,arg)


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
            fig,ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['ED']]
            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
            ax[i].set(xlim=(0, max_x), ylim=(0, 1.5))
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='ED m', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
            ax[i].set_yticks(ticks=[0,0.5,1,1.5])
            ax[i].set_yticklabels(labels=[0,0.5,1,1.5], fontsize=sets['yticksSIZE'])
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
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
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])

        fig.subplots_adjust(hspace=sets['hspace'], bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
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
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, 3))
    for i in range(len(chrom)):
        d=df[df['#CHROM'] == chrom[i]]
        max_x = int(arg['contigs'][chrom[i]])
        x=d[['POS']]
        y=d[['ED']]
        ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
        ax[i].set(xlim=(0, max_x), ylim=(0, 1.5))
        ax[i].set_xticks([])
        ax[i].set_xlabel(xlabel=f'{chrom[i]}', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        ax[i].set_ylabel(ylabel='ED m', fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'], rotation=90)
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

    fig.subplots_adjust(wspace=0, bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
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
            fig,ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
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
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
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
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)',  fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])

        fig.subplots_adjust(hspace=sets['hspace'], bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
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
    fig, ax = plt.subplots(1, len(chrom), figsize=(af, 3))
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
    fig.subplots_adjust(wspace=0, bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
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
    ylab = 'A. F.'
    boost = False
    multi = (True,False)    #Multivertical Not Manhattan
    labs_list = list()
    f_name = list()
    arg, ticks_y, lim_y, ylab, t, rt, key2, key = chooseOptionsAFPlots(arg, multi, g_type, ticks_y, lim_y, ylab)

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

            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
            ax[i].set(xlim=(0, max_x), ylim=lim_y)
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel=ylab, fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x = sets['titleXPOS'], y=sets['titleYPOS'])
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
                if arg['--boost'] != False:
                    keyB = 36
                    if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                        plot_boost(d, arg, ax[i], max_x)
                        boost = True
                    if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                        plot_boost(d, arg, ax[i], max_x)
                        boost = True
                if arg['--distance-boost'] != False:
                    keyB = 35
                    if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                        plotDistanceBoost(d, arg, ax[i], max_x, chrom[i])
                        boost = True
                    if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                        plotDistanceBoost(d, arg, ax[i], max_x, chrom[i])
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
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=sets['xlabSIZE'], labelpad=sets['xlabDIST'])
        fig.subplots_adjust(hspace=sets['hspace'], bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
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
                if arg['--boost'] != False:
                    cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--boost'])))
                if arg['--distance-boost'] != False:
                    cap.append(arg['lines'][keyB].format(arg['color_names']['BOOST'], str(arg['--distance-boost'])))
        if arg['--ci95'] and g_type == 'DELTA':
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
            fig, ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
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
            ax[i].set_ylabel(ylabel='A. F.', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=sets['titleSIZE'], rotation=0, x =sets['titleXPOS'], y=sets['titleYPOS']) 
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))  
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
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=sets['xlabSIZE'],labelpad=sets['xlabDIST'])

        fig.subplots_adjust(hspace=sets['hspace'], bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
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
            fig,ax=plt.subplots(len(chrom), 1, figsize=(10, 2.5*len(chrom)))
        if len(chrom) == 1:
            fig, ax=plt.subplots(2,1,figsize=(10, 2.5*2))
        for i in range(len(chrom)):
            d=df[df['#CHROM'] == chrom[i]]
            x=d[['POS']]
            y=d[['SNPidx2']]
            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
            ax[i].set(xlim=(0, max_x), ylim=lim_y)
            ax[i].set_xticks(ticks=np.arange(0, max_x, 5e6))
            ax[i].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
            ax[i].tick_params(labelbottom=False)
            ax[i].set_ylabel(ylabel='A. F.', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
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
                ax[i].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=sets['xlabSIZE'],labelpad=sets['xlabDIST'])

        fig.subplots_adjust(hspace=sets['hspace'], bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
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
        ax[0].scatter(x, y1, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        ax[1].scatter(x, y2, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        ax[2].scatter(x, y3, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
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
        if arg['--moving-avg'] != False:
            k1 = 3
            k2 = 5
        if arg['--distance-avg'] != False:
            k1 = 2
            k2 = 4
    else:
        key = 31
        if arg['--moving-avg'] != False:
            k1 = 9
            k2 = 11
        if arg['--distance-avg'] != False:
            k1 = 8
            k2 = 10
    for i in range(len(chrom)):
        fig, ax=plt.subplots(3, 2, figsize=(8.5, 8.5))#(7,9)paper
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
        ax[0,0].set_ylabel('SNP-index 1', fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'])
        ax[0,0].yaxis.set_label_coords(sets['ylabXPOS'],sets['ylabYPOS'])
        ax[0,0].set_title('(a)', fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
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
        ax[1,0].set_ylabel('SNP-index 2', fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'])
        ax[1,0].yaxis.set_label_coords(sets['ylabXPOS'],sets['ylabYPOS'])
        ax[1,0].set_title('(b)', fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
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
        ax[2,0].set_title('(c)', fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
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
        ax[2,0].set_xlabel(xlabel='Chromosomal position (bp)',fontsize=sets['xlabSIZE'],labelpad=sets['xlabDIST'])
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
        ax[0,1].set_ylabel(ylabel='Euclidean distance', fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'], rotation=90)
        ax[0,1].yaxis.set_label_coords(sets['ylabXPOS'],sets['ylabYPOS'])
        ax[0,1].set_title('(d)',fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
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
        ax[1,1].set_title('(e)', fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
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
        ax[2,1].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[2,1].set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=sets['xticksSIZE'])
        ax[2,1].tick_params(axis='x', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        ax[2,1].set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)',fontsize=sets['ylabSIZE'],labelpad=sets['ylabDIST'])
        ax[2,1].yaxis.set_label_coords(sets['ylabXPOS'],sets['ylabYPOS'])
        ax[2,1].set_title('(f)', fontsize=sets['titleSIZE'], rotation=0, x=sets['titleXPOS'], y=sets['titleYPOS'])
        ax[2,1].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[2,1], 'log10PVALUE')
        if arg['--distance-avg'] != False:
            plotDistanceAvg(d, chrom[i], arg, ax[2,1], 'log10PVALUE')
        if arg['--bonferroni']:
            max_x_ch=(max(d['POS'])/max_x)
            ax[2,1].axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        ax[2,1].xaxis.set_major_formatter(ScalarFormatter())
        ax[2,1].ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax[2,1].xaxis.get_offset_text().set_fontsize(sets['xSCIoffsetSIZE'])
        ax[2,1].set_xlabel(xlabel='Chromosomal position (bp)',fontsize=sets['xlabSIZE'],labelpad=sets['xlabDIST'])
        ax[2,1].spines['top'].set_visible(False)
        ax[2,1].spines['right'].set_visible(False)
        ax[2,1].tick_params(axis='y', which='major', labelsize=sets['yticksSIZE'])
        ax[2,1].tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
        for spine in ax[2,0].spines.values():
            spine.set_linewidth(sets['axSpinesWIDTH'])
        # Save file
        fig.subplots_adjust(hspace=sets['hspace'], wspace=sets['wspace'], bottom=sets['bottomDIST'],top=sets['topDIST'],right=sets['rightDIST'],left=sets['leftDIST'])
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
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][key].format(
                                        arg['lines'][k1].format(arg['color_names']['SNPidx1'], str(arg['--moving-avg'])),\
                                        arg['lines'][k2].format(arg['color_names']['SNPidx2'], str(arg['--moving-avg'])),\
                                        '' if arg['--ref-genotype'] == True else ', in absolute value',\
                                        arg['lines'][22].format(arg['color_names']['ci'], str(arg['n_markers'])) if arg['--ci95'] else '',\
                                        arg['lines'][24].format(arg['color_names']['mvg']) if isinstance(ed100,pd.DataFrame) else '',\
                                        arg['lines'][15].format(str(arg['n_markers'])) if arg['--bonferroni'] else ''
                                        ))
            
                cap.append(arg['lines'][29].format(arg['color_names']['mvg'], str(arg['--moving-avg'])))
            elif arg['--distance-avg'] != False:
                cap.append(arg['lines'][key].format(
                                        arg['lines'][k1].format(arg['color_names']['SNPidx1'], str(arg['--distance-avg'])),\
                                        arg['lines'][k2].format(arg['color_names']['SNPidx2'], str(arg['--distance-avg'])),\
                                        '' if arg['--ref-genotype'] == True else ', in absolute value',\
                                        arg['lines'][22].format(arg['color_names']['ci'], str(arg['n_markers'])) if arg['--ci95'] else '',\
                                        arg['lines'][24].format(arg['color_names']['mvg']) if isinstance(ed100,pd.DataFrame) else '',\
                                        arg['lines'][15].format(str(arg['n_markers'])) if arg['--bonferroni'] else ''
                                        ))
                cap.append(arg['lines'][28].format(arg['color_names']['mvg'], str(arg['--distance-avg'])))
            else:
                cap.append(arg['lines'][key].format(
                                        '',\
                                        '',\
                                        '' if arg['--ref-genotype'] == True else ', in absolute value',\
                                        arg['lines'][22].format(arg['color_names']['ci']) if arg['--ci95'] else '',\
                                        arg['lines'][24].format(arg['color_names']['mvg']) if isinstance(ed100,pd.DataFrame) else '',\
                                        arg['lines'][15].format(str(arg['n_markers'])) if arg['--bonferroni'] else ''
                                        ))
            write_caption(f,cap, arg)

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
        ax[0].scatter(x, y1, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        ax[1].scatter(x, y2, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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

        ax[2].scatter(x, y3, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        if arg['--ci95'] and (arg['--moving-avg'] != False or arg['--distance-avg'] != False):
            if arg['--moving-avg'] != False:
                calc_ci(d, arg, ax[2])
            if arg['--distance-avg'] != False:
                distanceCI(d, arg, ax[2], chrom[i])
            ax[2].axhline(y=0, color = 'black', linestyle='dashed', linewidth=0.75)
        ax[2].spines['top'].set_visible(False)
        ax[2].spines['right'].set_visible(False)

        # log10PVALUE
        ax[3].scatter(x, y4, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
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
        ax[0].scatter(x, y1, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        ax[1].scatter(x, y2, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        ax[2].scatter(x, y3, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
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
    ax[0].plot(d['mediamovilx'], d['mediamovily1'], c=arg['--palette']['SNPidx1'], lw=arg['LINE_W'])
    #SNPidx2
    d['producty2']=d['SNPidx2'] * d['DP2']
    d['mediamovily2']=(d['producty2'].rolling(arg['--moving-avg']).sum() / d['DP2'].rolling(arg['--moving-avg']).sum())
    ax[1].plot(d['mediamovilx'], d['mediamovily2'], c=arg['--palette']['SNPidx2'], lw=arg['LINE_W'])
    # D-SNPidx
    d['productyD']=d['DELTA'] * d['DP']
    d['mediamovilD']=(d['productyD'].rolling(arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())
    ax[2].plot(d['mediamovilx'], d['mediamovilD'], c=arg['--palette']['DELTA'], lw=arg['LINE_W'])


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
    ax[0].plot(res['POSx'], res['avg1'], c=arg['--palette']['SNPidx1'], lw=arg['LINE_W'])
    ax[1].plot(res['POSx'], res['avg2'], c=arg['--palette']['SNPidx2'], lw=arg['LINE_W'])
    ax[2].plot(res['POSx'], res['avgD'], c=arg['--palette']['DELTA'], lw=arg['LINE_W'])


def plot_boostManhattan(d, arg, ax, max_x, chrom, i):
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
    c_=arg['--palette']['BOOST']
    ax2=ax.twinx()
    ax2.spines['top'].set_visible(False)
    ax2.set(xlim=(0, max_x), ylim=(0, 1))
    ax2.set_yticks([])
    ax2.plot(d['medboostx'], d['medboost'], c=c_, lw=arg['LINE_W'], linestyle='dashed')
    

    if chrom[i] == chrom[-1]:
        ax2.axes.get_yaxis().set_visible(True)
        ax2.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
        ax2.set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=man_yticksS)
        ax2.set_ylabel(ylabel='Boost', fontsize=man_ylabS, rotation=90, labelpad=man_ylabD)

def plot_boost(d, arg, ax, max_x, sets):
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
    ax2.plot(d['medboostx'], d['medboost'], c=c_, lw=arg['LINE_W'], linestyle='dashed')
    ax2.set(xlim=(0, max_x), ylim=(0, 1))
    ax2.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
    ax2.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
    ax2.set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=sets['yticksSIZE'])
    ax2.set_ylabel(ylabel='Boost', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])

def plotDistanceBoostManhattan(d, arg, ax, max_x, chrom, i):
    # Mediamovil Boost
    if 'SNPidx2' not in arg['--fields']:
        d['DP']=d['DPdom_1']+d['DPrec_1']
    else:
        d['DP']=d['DPdom_2']+d['DPrec_2']
    d['BOOST']=d['BOOST'] * arg['lim']
    d['prodx'] = d['DP']*d['POS']
    d['prody'] = d['DP']*d['BOOST']
    interval_ranges = [x for x in range(0, max_x, arg['--distance-boost'])]
    if interval_ranges[-1] < max_x:
        interval_ranges.append(max_x)
    interval_labels = [f'{start}-{end-1}' for start, end in zip(interval_ranges, interval_ranges[1:])]
    d['Interval'] = pd.cut(d['POS'], bins=interval_ranges, labels=interval_labels)
    res = d.groupby('Interval').agg({'prodx': sum, 'prody': sum, 'DP': sum})
    res['POSx'] = (res['prodx']/res['DP'])
    res['VALy'] = (res['prody']/res['DP'])
    res = res.dropna()

    c_=arg['--palette']['BOOST']
    ax2=ax.twinx()
    ax2.spines['top'].set_visible(False)
    ax2.set(xlim=(0, max_x), ylim=(0, 1))
    ax2.set_yticks([])
    ax2.plot(res['POSx'], res['VALy'], c=c_, lw=arg['LINE_W'], linestyle='dashed')
    
    
    if chrom[i] == chrom[-1]:
        ax2.axes.get_yaxis().set_visible(True)
        ax2.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
        ax2.set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=man_yticksS)
        ax2.set_ylabel(ylabel='Boost', fontsize=man_ylabS, rotation=90, labelpad=man_ylabD)


def plotDistanceBoost(d, arg, ax, max_x, chrom, sets):
    # Mediamovil Boost
    if 'SNPidx2' not in arg['--fields']:
        d['DP']=d['DPdom_1']+d['DPrec_1']
    else:
        d['DP']=d['DPdom_2']+d['DPrec_2']
    d['BOOST']=d['BOOST'] * arg['lim']
    max_x = arg['contigs'][chrom]
    d['prodx'] = d['DP']*d['POS']
    d['prody'] = d['DP']*d['BOOST']
    interval_ranges = [x for x in range(0, max_x, arg['--distance-boost'])]
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
    ax2.plot(res['POSx'], res['VALy'], c=c_, lw=arg['LINE_W'], linestyle='dashed')
    ax2.set(xlim=(0, max_x), ylim=(0, 1))
    ax2.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
    ax2.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
    ax2.set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=sets['yticksSIZE'])
    ax2.set_ylabel(ylabel='Boost', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])

def plot_ED100_4(arg, ax, max_x, max_y, ed100, ch, sets):
    flag = ch in ed100['#CHROM'].values
    ax2 = ax.twinx()
    c_ = arg['--palette']['mvg']
    ax2.set(xlim=(0, max_x), ylim=(0, max_y))
    ax2.set_yticks(ticks=[0, 0.5, 1, 1.5])
    ax2.set_yticklabels(labels=[0, 0.5, 1, 1.5], fontsize=sets['yticksSIZE'])
    ax2.tick_params(axis='y', which='both', length=sets['ticksSpinesLENGTH'], width=sets['ticksSpinesWIDTH'])
    #ax2.tick_params(axis='y', which='major', labelsize=10)
    if flag:
        ed100ch = ed100[ed100['#CHROM'] == ch]
        #print(ed100ch)
        ax2.set_ylabel(ylabel='ED $\mathregular{100^4}$  $\mathregular{x10^8}$ ', fontsize=sets['ylabSIZE'], rotation=90, labelpad=sets['ylabDIST'])
        ax2.plot(ed100ch['POS'], ed100ch['ED100_4'], c=c_, lw=arg['LINE_W'])
    ax2.spines['top'].set_visible(False)

def create_caption(arg, res_tit):
    f= io.open(arg['captions_dir'] + res_tit + '.txt', 'a', encoding='utf8')
    return f

def create_subfolder(arg, sf):
    wd = os.getcwd()
    arg['sub_folder'] = wd+'/'+arg['--outdir']+'/'+sf
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
            ax[i].scatter(x, y, s=arg['DOT_SIZE'], c=arg['--palette']['dots'], alpha=arg['--alpha'], clip_on=False)
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
        plt.savefig(arg['sub_folder']+filename, dpi=arg['DPI'], transparent=True)
        plt.close()