from pylibs.analysis_lib import *
from pylibs.constants import *
import warnings
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging
import scipy.stats as st
import sys
import os
from math import log
import errno
from docopt import docopt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter,
                               AutoMinorLocator)

# warnings.filterwarnings("ignore")
# sys.path.append(os.getcwd())
pd.options.mode.chained_assignment = None
# ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','MAX_SNPidx2','BOOST']
# ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','SNPidx1','BOOST']
# ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','SNPidx1','SNPidx2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE']
# ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','MAX_SNPidx2','FISHER','BOOST','PVALUE','log10PVALUE']
# ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','PVALUE','log10PVALUE']
# ['#CHROM','POS','REF','ALT','DPref_1','DPalt_1','DPref_2','DPalt_2','SNPidx1','SNPidx2','DELTA','PVALUE','log10PVALUE']
fields = {'#CHROM': str, 'POS': int, 'REF': str, 'ALT': str, 'DPref_1': int, 'DPalt_1': int, 'DPref_2': int, 'DPalt_2': int, 'SNPidx1': float,
          'SNPidx2': float, 'DELTA': float, 'MAX_SNPidx2': float, 'PVALUE': float, 'log10PVALUE': float, 'FISHER': float, 'BOOST': float, 'ED': float, 'G':float, 'CI95': float}


def test_plot(arg):
    font_files = matplotlib.font_manager.findSystemFonts()
    for font_file in font_files:
        matplotlib.font_manager.fontManager.addfont(font_file)
    plt.rcParams['font.family'] = 'Arial'
    global fields
    global header
    inp_f = arg['<input_file>']
    inp_ext = inp_f.split('.')[-1]
    if inp_ext == 'csv':
        sep_ = ','
    else:
        sep_ = '\t'
    try:
        df = pd.read_csv(inp_f, sep=sep_, dtype=fields, header='infer')
    except FileNotFoundError:
        print('Error: The input file {} does not exist'.format(inp_f))
        sys.exit()
    wd = os.getcwd()
    try:
        os.makedirs(wd+'/'+arg['--outdir'])
    except FileExistsError:
        print('Warning: the output directory already exists')
        pass
    arg['--outdir'] = wd+'/'+arg['--outdir']+'/'
    if arg['--captions']:
        arg['captions_dir'] = arg['--outdir']+'captions/'
        try:
            os.makedirs(arg['captions_dir'])
        except FileExistsError:
            pass

    chroms = df['#CHROM'].unique()
    if arg['--chromosomes'] == 'all':
        arg['--chromosomes'] = list(chroms)
    else:
        ch_idx = [int(c) - 1 for c in arg['--chromosomes'].split(',')]
        chroms_ = [chroms[idx] for idx in ch_idx]
        arg['--chromosomes'] = chroms_

    arg['labs'] = {arg['--chromosomes'][i]:chr(97+i) for i in range(len(arg['--chromosomes']))}
    chrom_lists = [arg['--chromosomes'][i:i+6] for i in range(0, len(arg['--chromosomes']), 6)]
    arg['chrom_lists'] = chrom_lists
    if arg['--multi-chrom'] == True and len(arg['--chromosomes']) == 1:
        arg['--multi-chrom'] = False
    arg['--fields'] = list(df.columns)
    arg['n_markers'] = len(df)
    arg['--alpha'] = float(arg['--alpha'])
    arg['lim'] = 1e-90
    arg['--fileformat'] = '.'+arg['--fileformat']
    if arg['--moving-avg'] == None:
        arg['--moving-avg'] = False
    else:
        arg['--moving-avg'] = int(arg['--moving-avg'])
        if arg['--moving-avg'] <= 0:
            print(
                'Error: The window size  for m. avg must be an integer higher than zero.')
            sys.exit()
    if 'log10PVALUE' in arg['--fields']:
        r = pd.DataFrame(dtype=np.float32)
        r['log10PVALUE'] = df['log10PVALUE']
        r = r.replace({-np.inf: r[np.isfinite(r)].min().min()})
        # r = r.replace({np.nan: r[np.isfinite(r)].min().min()})
        # print(r[np.isfinite(r)].min().min())
        df['log10PVALUE'] = r['log10PVALUE']
    # Checking graphic types
    if 'mbsplot' in arg.keys():
        arg = check_mbs_opts(arg)
    if 'qtlplot' in arg.keys():
        arg = check_qtl_opts(arg)
    return df, arg


def test_merge(arg):
    global fields
    inp_f = arg['<input_file>']
    inp_ext = inp_f.split('.')[-1]
    if inp_ext == 'csv':
        sep_ = ','
    else:
        sep_ = '\t'
    try:
        df = pd.read_csv(inp_f, sep=sep_, dtype=fields, header='infer')
    except FileNotFoundError:
        print('Error: The input file {} does not exist'.format(inp_f))
        sys.exit()
    wd = os.getcwd()
    try:
        os.makedirs(wd+'/'+arg['--outdir'])
    except FileExistsError:
        print('Warning: the output directory already exists')
        pass

    arg['--outdir'] = wd+'/'+arg['--outdir']+'/'
    arg['--window'] = int(arg['--window'])
    if arg['--window'] <= 0:
        print('Error: The window size must be an integer higher than zero.')
        sys.exit()
    if arg['--output'] != None:
        arg['--fileformat'] = '.'+arg['--fileformat']
        arg['--output'] = str(arg['--window'])+'w_'+arg['--output'] + arg['--fileformat']
        arg['--output'] = check_save_an(arg, arg['--output'])
    else:
        arg['--output'] = None
    if arg['--fileformat'] == '.csv':
        arg['spacer'] = ','
    else:
        arg['spacer'] = '\t'

    chroms = df['#CHROM'].unique()
    if arg['--chromosomes'] == 'all':
        arg['--chromosomes'] = list(chroms)
    else:
        ch_idx = [int(c) - 1 for c in arg['--chromosomes'].split(',')]
        chroms_ = [chroms[idx] for idx in ch_idx]
        arg['--chromosomes'] = chroms_
    arg['--fields'] = list(df.columns)
    arg['--fields'].remove('REF')
    arg['--fields'].remove('ALT')
    arg['--fileformat'] = '.'+arg['--fileformat']
    arg['--window'] = int(arg['--window'])
    if arg['--window'] <= 0:
        print('Error: The window size must be an integer higher than zero.')
        sys.exit()
    # Checking graphic types
    return df, arg


def check_mbs_opts(arg):
    if arg['--boost'] == None:
        arg['--boost'] = False
    else:
        arg['--boost'] = int(arg['--boost'])
        if arg['--boost'] <= 0:
            print('Error: The window size for boost must be an integer higher than zero.')
            sys.exit()
    arg['titles'] = titles_mbs
    arg['lines'] = lines_mbs
    palette = palettes_mbs[arg['--palette']]
    arg['--palette'] = {k: v[0] for k, v in palette.items()}
    arg['color_names'] = {k: v[1] for k, v in palette.items()}
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
        if 'ED' in arg['--fields']:
            arg['--euclidean-distance'] = True
        if 'MAX_SNPidx2' in arg['--fields'] and 'SNPidx1' not in arg['--fields'] and 'SNPindx2' not in arg['--fields']:
            arg['--max-allele-freq2'] = True
        else:
            arg['--max-allele-freq2'] = False
    else:
        if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields'] and arg['--combine'] == True:
            print('Error: is not possible make combined graphics with your input data')
            sys.exit()
        if arg['--pvalue'] == True and 'log10PVALUE' not in arg['--fields']:
            print('Error: is not possible make P-VALUE graphics with your input data')
            sys.exit()
        if arg['--euclidean-distance'] == True and 'ED' not in arg['--fields']:
            print('Error: is not possible make Euclidean Distance graphics with your input data')
            sys.exit()
        if arg['--allele-freq-1'] == True and 'SNPidx1' not in arg['--fields']:
            print('Error: is not possible make phased Allele Frequency graphics with your input data. Please use -M or -a option')
            sys.exit()
        if arg['--allele-freq-2'] == True and 'SNPidx2' not in arg['--fields']:
            print('Error: is not possible make Allele Frequency graphics for this pool. Please use -R, -M or -a option')
            sys.exit()
        if arg['--max-allele-freq2'] == True and 'MAX_SNPidx2' not in arg['--fields']:
            print('Error: is not possible make phased MAX Allele Frequency graphics with your input data. Please use -D or -R or -a option')
            sys.exit()
    if arg['--all'] == False and arg['--pvalue'] == False and arg['--allele-freq-1'] == False and arg['--allele-freq-2'] == False:
        print('Warning: any graphic was selected. Please use -a option to make all possible graphics with your data.')
    return arg


def check_qtl_opts(arg):
    arg['titles'] = titles_qtl
    arg['lines'] = lines_qtl
    palette = palettes_qtl[arg['--palette']]
    arg['--palette'] = {k: v[0] for k, v in palette.items()}
    arg['color_names'] = {k: v[1] for k, v in palette.items()}
    if arg['--all']:
        if 'log10PVALUE' in arg['--fields']:
            arg['--pvalue'] = True
        else:
            arg['--pvalue'] = False
        if 'SNPidx1' in arg['--fields'] and 'SNPidx2' in arg['--fields']:
            arg['--allele-freq-H'] = True
            arg['--allele-freq-L'] = True
            arg['--combine'] = True
        if 'DELTA' in arg['--fields']:
            arg['--delta'] = True
            arg['--ci95'] = True
        if len(arg['--chromosomes']) > 1:
            arg['--multi-chrom'] = True

    else:
        if arg['--pvalue'] == True and 'log10PVALUE' not in arg['--fields']:
            print('Error: is not possible make P-VALUE graphics with your input data')
            sys.exit()
        if arg['--delta'] == True and 'DELTA' not in arg['--fields']:
            print(
                'Error: is not possible make phased SNP-idx graphics with your input data. Please use -a.')
            sys.exit()
        if arg['--ci95'] == True and 'DELTA' not in arg['--fields']:
            print('Error: is not possible to calculate confidense interval without DELTA field. Please use -a.')
            sys.exit()
        if ('SNPidx1' not in arg['--fields'] or 'SNPidx2' not in arg['--fields']) and (arg['--allele-freq-H'] == True or arg['--allele-freq-L'] == True):
            print('Error: is not possible make phased SNP-idx graphics with your input data. Please use -a or -p.')
            sys.exit()
        if arg['--combine'] == True and 'SNPidx1' not in arg['--fields'] or 'SNPidx2' not in arg['--fields']:
            print('Error: is not possible to make combined graphics with your input data. Please use -a.')
            sys.exit()
    if arg['--all'] == False and arg['--pvalue'] == False and arg['--delta'] == False and arg['--allele-freq-H'] == False and arg['--allele-freq-L'] == False:
        print('Warning: any graphic was selected. Please use -a option to make all possible graphics with your data.')

    return arg


def check_save(arg, file_name):
    typ = arg['--fileformat']
    if os.path.isfile(arg['--outdir']+file_name):

        expand = 0
        while True:
            expand += 1
            nw_file_name = file_name.split(typ)[0] + '(' + str(expand) + ')' + typ
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
        d.index = np.arange(len(d))
        d_list = [i for i in range(len(d))]
        g_list = [d_list[i:i+w]
            for i in range(0, len(d_list), w)]  # indexes for each group

        for g in g_list:  # for each group in group list
            res = list()
            if 'DPref_2' not in arg['--fields'] and 'DPalt_2' not in arg['--fields']:
                rAt, rBt, rPost, rTt = 0, 0, 0, 0
            else:
                rAt, rBt, rCt, rDt, rPost, rTt = 0, 0, 0, 0, 0, 0
            for i in g:  # for each index in group
                if 'DPref_2' not in arg['--fields'] and 'DPalt_2' not in arg['--fields']:
                    rA = d.loc[i].DPref_1
                    rB = d.loc[i].DPalt_1
                    rT = rA + rB
                    rPos = d.loc[i].POS*rT
                    rTt += rT
                    rAt += rA
                    rBt += rB
                    rPost += rPos
                else:
                    rA = d.loc[i].DPref_1
                    rB = d.loc[i].DPalt_1
                    rC = d.loc[i].DPref_2
                    rD = d.loc[i].DPalt_2
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
            if 'DPref_2' in arg['--fields'] and 'DPalt_2' in arg['--fields']:
                res += [rCt, rDt]
                if 'SNPidx1' in arg['--fields'] and 'SNPidx2' in arg['--fields']:
                    rSNPidx1 = rBt/(rAt+rBt)
                    rSNPidx2 = rDt/(rCt+rDt)
                    res += [rSNPidx1, rSNPidx2]
                if 'MAX_SNPidx2' in arg['--fields']:
                    rMAX_SNPidx2 = (max(rDt, rCt))/(rDt+rCt)
                    fisher = LogFisher(rAt, rBt, rCt, rDt)
                    boost = 1/(sys.float_info.min +
                               abs(1 - 1/max(rSNPidx2, 1-rSNPidx2)))
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
    chrom = arg['--chromosomes']
    ED100 = np.array([])
    for ch in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[ch]]
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
    df['ED100_4'] = pd.Series(ED100**4)
    return df


def plot_ED(df, arg):
    rang = 100
    df = get_ED100_4(df, arg, rang)
    chrom = arg['--chromosomes']
    typ = arg['--fileformat']
    max_y = max(df['ED100_4'])
    for ch in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[ch]]
        d.index = np.arange(len(d))
        max_x = max(d['POS'])
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
        ax.set_ylabel(ylabel='EDm', fontsize=12, rotation=90, labelpad=15)
        ax.set_yticks(ticks=[0, 0.5, 1, 1.5])
        ax.set_yticklabels(labels=[0, 0.5, 1, 1.5], fontsize=8)
        ax.spines['top'].set_visible(False)
        if arg['--moving-avg'] != False:
            ax2 = ax.twinx()
            c_ = arg['--palette']['mvg']
            ax2.set(xlim=(0, max_x), ylim=(0, max_y))
            ax2.tick_params(axis='y', which='major', labelsize=8)
            ax2.yaxis.set_major_formatter(ScalarFormatter())
            ax2.ticklabel_format(axis='y', style='scientific', scilimits=(7,8), useMathText=True)
            ax2.set_ylabel(ylabel='ED$\mathregular{100^4}$', fontsize=12, rotation=90, labelpad=15)
            ax2.plot(d['POS'], d['ED100_4'], c=c_, lw=1.25)
            ax2.spines['top'].set_visible(False)

        filename = 'ED_chr{}_{}'.format(chrom[ch], 'ED100_4' if arg['--moving-avg'] != False else '') + typ
        filename = check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename)
        plt.close()

def plot_G(df, arg):
    chrom = arg['--chromosomes']
    typ = arg['--fileformat']
    min_y, max_y =min(df['G']), max(df['G'])
    for ch in range(len(chrom)):
        d = df[df['#CHROM'] == chrom[ch]]
        max_x = max(d['POS'])
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

        filename = 'G_chr{}'.format(chrom[ch]) + typ
        filename = check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename)
        plt.close()

def pval_multi_graph(df, arg):
    print ("p_valor_multigrafico")
    t, rt = arg['titles'][1]
    if arg['--captions']:
        f = create_caption(arg, rt)
    typ = arg['--fileformat']
    min_y = min(df['log10PVALUE'])*1.05
    labs_list = list()
    f_name = list()
    for chrom_list in arg['chrom_lists']:
        chrom = chrom_list
        if len(chrom) > 1:
            fig, ax = plt.subplots(1, len(chrom), figsize=(len(chrom)*(3.3333333335), 2))
        if len(chrom) == 1:
            fig, ax = plt.subplots(1,2,figsize=(2*(3.3333333335), 2))
        for i in range(len(chrom)):
            d = df[df['#CHROM'] == chrom[i]]
            max_x = max(d['POS'])
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
        plt.savefig(arg['--outdir']+filename)
        plt.close()
    cap = list()
    if arg['--captions']:
        cap.append(', '.join(f_name))
        cap.append(t.format(', '.join(labs_list)))
        cap.append(arg['lines'][0].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            cap.append(arg['lines'][1].format(arg['color_names']['log10PVALUE'], str(arg['--moving-avg'])))
        if arg['--bonferroni']:
            cap.append(arg['lines'][2].format(str(arg['n_markers'])))
        write_caption(f, cap)
    
def pval_mono_graph(df, arg):
    chrom=arg['--chromosomes']
    typ=arg['--fileformat']
    
    min_y=min(df['log10PVALUE'])*1.05
    
    for i in range(len(chrom)):
        d=df[df['#CHROM'] == chrom[i]]
        max_x=max(d['POS'])
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

        plt.savefig(arg['--outdir']+filename)
        plt.close()

        if arg['--captions']:
            f = create_caption(arg, rtch)
            cap.append(filename)
            cap.append(t.format(chrom[i]))
            cap.append(arg['lines'][0].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][1].format(arg['color_names']['log10PVALUE'], str(arg['--moving-avg'])))
            if arg['--bonferroni']:
                cap.append(arg['lines'][2].format(str(arg['n_markers'])))
            write_caption(f, cap)


def qq_plot(df, arg):
    #TODO
    df['DP']=df['DPref_1']+df['DPalt_1']+df['DPref_2']+df['DPalt_2']
    chrom=arg['--chromosomes']
    typ=arg['--fileformat']
    inp_file=arg['<input_file>']
    r=pd.DataFrame()
    r['A']=df.DPref_1
    r['B']=df.DPalt_1
    r['C']=df.DPref_2
    r['D']=df.DPalt_2
    r['#CHROM']=df['#CHROM']
    r['obs_logPVALUE']=-1*df.log10PVALUE
    r['e_logPVALUE']=r.apply(expected_pvalue, axis=1)
    r=r.dropna()
    for i in range(len(chrom)):
        r_=r[r['#CHROM'] == chrom[i]]
        x=r_['e_logPVALUE']
        y=r_['obs_logPVALUE']
        fig, ax=plt.subplots(figsize=(5, 5))
        ax.scatter(x, y, s=0.5, c='royalblue', vmin=0, vmax=100, alpha=0.4)
        # ax.set(xlim=(0, max_x), ylim=(min_y, 0))
        # ax.set_xticks(ticks=np.arange(0,max_x,5e6))#,fontdict={'family':'arial', 'size':'8'})
        # ax.set_xticklabels(labels=np.arange(0,max_x,5e6),fontdict={'family':'arial', 'size':'8'})
        # ax.xaxis.set_major_formatter(ScalarFormatter())
        # ax.ticklabel_format(axis='x', style='scientific', scilimits=(6,6), useMathText=True)
        ax.set_xlabel(xlabel='Expected '+'-log' +
                      r'$_{10}$'+'(p-value)', fontdict={'family': 'arial', 'size': '8'})
        ax.set_ylabel(ylabel='Observed '+'log' +
                      r'$_{10}$'+'(p-value)', fontdict={'family': 'arial', 'size': '8'})
        # ax.set_yticks(ticks=np.arange(round(min_y),0,2))
        # ax.set_yticklabels(labels=np.arange(round(min_y),0,2),fontdict={'family':'arial', 'size':'8'})
        # ax.axhline(y=threshold, color='black', xmin=0,xmax=max_x_ch, linestyle='dashed', linewidth=0.75)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # mov_avg.to_csv(path_or_buf='./mov_avg_pvalue.txt',sep='\t', header=True)
        filename='chr_{}_qqplot'.format(chrom[i]) + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename)


def AF1_AF2_mono_graph(df, arg):
    chrom=arg['--chromosomes']
    typ=arg['--fileformat']

    for i in range(len(chrom)):
        d=df[df['#CHROM'] == chrom[i]]
        t1,rt1,_=arg['titles'][3]
        t2,rt2,_=arg['titles'][4]
        max_x=max(d['POS'])
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
            plot_avg(d, arg, ax, 'SNPidx1')
            plot_avg(d, arg, ax, 'SNPidx2')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        rtch1 = rt1.format(chrom[i])
        filename=rtch1 + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename)
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
            write_caption(f,cap) 
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
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        rtch2 = rt2.format(chrom[i])
        filename=rtch2 + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename)
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
            write_caption(f,cap)


def AF_mono_graph(df, arg, g_type):
    chrom=arg['--chromosomes']
    typ=arg['--fileformat']
    ticks_y=[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    ylab = 'Allele Frequency'
    boost = False
    for i in range(len(chrom)):
        if g_type == 'SNPidx1':
            key = 3
            key2 = -2
            t,_,rt=arg['titles'][key]
        if g_type == 'SNPidx2':
            key = 4
            key2 = -1
            t,_,rt=arg['titles'][key]
        if g_type == 'MAX_SNPidx2':
            ylab = 'Maximum Allele Frequency'
            key = 5
            key2 = -3
            t,rt=arg['titles'][key]
            ticks_y=[0.5, 0.75, 1]
            lim_y=(0.5, 1)
        if g_type == 'DELTA':
            ylab = '$\Delta$'+' (SNP-index)'
            key = 5
            key2 = -3
            t,rt = arg['titles'][key]
            ticks_y=[-1, 0, 1]
            lim_y=(-1.5, 1.5)
        d=df[df['#CHROM'] == chrom[i]]
        max_x=max(d['POS'])
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
        if 'mbsplot' in arg.keys():
            if arg['--boost'] != False:
                if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                    plot_boost(d, arg, ax, max_x)
                    boost = True
                if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                    plot_boost(d, arg, ax, max_x)
                    boost = True
        if 'qtlplot' in arg.keys():
            if arg['--ci95'] and arg['--moving-avg'] != False and g_type == 'DELTA':
                calc_ci(d, arg, ax)
                ax.axhline(y=0, color='black', linestyle='dashed', linewidth=0.75)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        rtch = rt.format(chrom[i])
        filename=rtch + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename)
        plt.close()
        cap = list()
        if arg['--captions']:
            f = create_caption(arg,rtch)
            cap.append(filename)
            cap.append(t.format(chrom[i]))
            cap.append(arg['lines'][key2].format(arg['color_names']['dots']))
            if arg['--moving-avg'] != False:
                cap.append(arg['lines'][key].format(arg['color_names'][g_type], str(arg['--moving-avg'])))
            if 'mbsplot' in arg.keys():
                if boost == True:
                    cap.append(arg['lines'][6].format(arg['color_names']['BOOST'], str(arg['--boost'])))
            if 'qtlplot' in arg.keys():
                if arg['--ci95'] and g_type == 'DELTA':
                    cap.append(arg['lines'][7].format(arg['color_names']['ci'], str(arg['n_markers'])))
            write_caption(f, cap)

def pval_multi_Vertical_graph(df, arg):
    typ=arg['--fileformat']
    min_y = min(df['log10PVALUE'])*1.05
    max_x = max(df['POS'])
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
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=18, rotation=0, x=-0.14, y=0.85)
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].tick_params(axis='y', which='major', labelsize=8)
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], 'log10PVALUE')
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
        plt.savefig(arg['--outdir']+filename)
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
        if arg['--bonferroni']:
            cap.append(arg['lines'][2].format(str(arg['n_markers'])))
        write_caption(f,cap)



def AF_multi_Vertical_graph(df, arg, g_type):
    typ=arg['--fileformat']
    max_x=max(df['POS'])
    ticks_y=[0, 0.25, 0.5, 0.75, 1]
    lim_y=(0, 1)
    ylab = 'Allele Frequency'
    boost = False
    labs_list = list()
    f_name = list()
    if g_type == 'SNPidx1':
        key = 7
        key2 = -2
        t,_,rt=arg['titles'][key]
    if g_type == 'SNPidx2':
        key = 8
        key2 = -1
        t,_,rt=arg['titles'][key]
    if g_type == 'MAX_SNPidx2':
        key = 9
        key2 = -3
        t,rt=arg['titles'][key]
        ylab='MAX Allele Frequency'
        ticks_y=[0.5, 0.75, 1]
        lim_y=(0.5, 1)
    if g_type == 'DELTA':
        ylab = '$\Delta$' + '(SNP-index)'
        key = 9
        key2 = -3
        t,rt = arg['titles'][key]
        ticks_y = [-1,0,1]
        lim_y = (-1.5, 1.5)
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
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=18, rotation=0, x = -0.14, y=0.85)
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))
            ax[i].set_yticks(ticks=ticks_y)
            ax[i].set_yticklabels(labels=ticks_y, fontsize=8)
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], g_type)
            if 'mbsplot' in arg.keys():
                if arg['--boost'] != False:
                    if 'SNPidx1' in arg['--fields'] and 'SNPidx2' not in arg['--fields']:
                        plot_boost(d, arg, ax[i], max_x)
                        boost = True
                    if g_type == 'SNPidx2' or g_type == 'MAX_SNPidx2':
                        plot_boost(d, arg, ax[i], max_x)
                        boost = True
            if 'qtlplot' in arg.keys():
                if arg['--ci95'] and arg['--moving-avg'] != False and g_type == 'DELTA':
                    calc_ci(d, arg, ax[i])
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
        plt.savefig(arg['--outdir']+filename)
        plt.close()
    cap = list()
    if arg['--captions']:
        f = create_caption(arg,rt)
        cap.append(', '.join(f_name))
        cap.append(t)
        cap = cap + labs_list
        cap.append(arg['lines'][key2].format(arg['color_names']['dots']))
        if arg['--moving-avg'] != False:
            key = key - 4
            cap.append(arg['lines'][key].format(arg['color_names'][g_type], str(arg['--moving-avg'])))
        if 'mbsplot' in arg.keys():
            if boost == True:
                cap.append(arg['lines'][6].format(arg['color_names']['BOOST'], str(arg['--boost'])))
        if 'qtlplot' in arg.keys() and g_type == 'DELTA':
            if arg['--ci95']:
                cap.append(arg['lines'][7].format(arg['color_names']['ci'], str(arg['n_markers'])))
        write_caption(f, cap)

def AF12_multi_Vertical_graph(df, arg):
    typ=arg['--fileformat']
    max_x=max(df['POS'])
    
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
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=18, rotation=0, x = -0.14, y=0.85) 
            labs_list.append('({}) Chromosome {}.'.format(arg['labs'][chrom[i]],chrom[i]))  
            ax[i].set_yticks(ticks=ticks_y)
            ax[i].set_yticklabels(labels=ticks_y, fontsize=8)
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], 'SNPidx2')
                plot_avg(d, arg, ax[i], 'SNPidx1')
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
        plt.savefig(arg['--outdir']+filename)
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
        write_caption(f,cap)
    
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
            ax[i].set_title('({})'.format(arg['labs'][chrom[i]]), fontsize=18, rotation=0, x = -0.14, y=0.85)
            ax[i].set_yticks(ticks=ticks_y)
            ax[i].set_yticklabels(labels=ticks_y, fontsize=8)
            if arg['--moving-avg'] != False:
                plot_avg(d, arg, ax[i], 'SNPidx1')
                plot_avg(d, arg, ax[i], 'SNPidx2')
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
        plt.savefig(arg['--outdir']+filename)
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
        write_caption(f,cap)

def AF1_AF2_pval_mono(df, arg):
    global fields
    chrom=arg['--chromosomes']
    typ=arg['--fileformat']
    min_y=min(df['log10PVALUE'])*1.05

    for i in range(len(chrom)):
        fig, ax=plt.subplots(3, 1, figsize=(8.5, 7))#(7,9)paper
        d=df[df['#CHROM'] == chrom[i]]
        max_x = max(d['POS'])
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
        ax[0].set_title('(a)', fontsize=17, rotation=0, x = -0.14, y=0.85)
        ax[0].set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
        ax[0].set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=8)
        if arg['--moving-avg'] != False:
            if arg['--combine']:
                plot_avg(d, arg, ax[0], 'SNPidx2')
                plot_avg(d, arg, ax[0], 'SNPidx1')
            else:
               plot_avg(d, arg, ax[0], 'SNPidx1') 
        ax[0].spines['top'].set_visible(False)
        ax[0].spines['right'].set_visible(False)

        # SNPidx2
        ax[1].scatter(x, y2, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[1].set(xlim=(0, max_x), ylim=(0, 1))
        ax[1].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[1].tick_params(labelbottom=False)
        ax[1].set_ylabel(ylabel='Allele Frequency', fontsize=12, rotation=90, labelpad=15)
        ax[1].set_title('(b)', fontsize=17, rotation=0, x = -0.14, y=0.85)
        ax[1].set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1])
        ax[1].set_yticklabels(labels=[0, 0.25, 0.5, 0.75, 1], fontsize=8)
        if arg['--moving-avg'] != False:
            if arg['--combine']:
                plot_avg(d, arg, ax[1], 'SNPidx1')
                plot_avg(d, arg, ax[1], 'SNPidx2')
            else:
               plot_avg(d, arg, ax[1], 'SNPidx2')
        ax[1].spines['top'].set_visible(False)
        ax[1].spines['right'].set_visible(False)

        # log10PVALUE
        ax[2].scatter(x, y3, s=0.5, c=arg['--palette']['dots'], alpha=arg['--alpha'])
        ax[2].set(xlim=(0, max_x), ylim=(min_y, 0))
        ax[2].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[2].set_xticklabels(labels=np.arange(0, max_x, 5e6),fontsize=8)
        ax[2].set_ylabel(ylabel='log'+r'$_{10}$'+'(p-value)',fontsize=12, labelpad=15)
        ax[2].set_title('(c)', fontsize=17, rotation=0, x = -0.14, y=0.85)
        ax[2].tick_params(axis='y', which='major', labelsize=8)
        if arg['--moving-avg'] != False:
            plot_avg(d, arg, ax[2], 'log10PVALUE')
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
        plt.savefig(arg['--outdir']+filename)
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
            if arg['--bonferroni']:
                cap.append(arg['lines'][2].format(str(arg['n_markers'])))
            write_caption(f,cap)


def snp_index_graph(df, arg):
    global fields
    chrom=arg['--chromosomes']
    typ=arg['--fileformat']
    t,rt=arg['titles'][11]

    for i in range(len(chrom)):
        fig, ax=plt.subplots(3, 1, figsize=(7.5, 5))
        d=df[df['#CHROM'] == chrom[i]]
        max_x=max(d['POS'])
        x=d[['POS']]
        y1=d['SNPidx1']
        y2=d['SNPidx2']
        y3=y2-y1  # delta

        #SNPidx1
        ax[0].scatter(x, y1, s=0.5, c=arg['--palette']['SNPidx1'], alpha=arg['--alpha'])
        ax[0].set(xlim=(0, max_x), ylim=(0, 1))
        ax[0].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[0].tick_params(labelbottom=False)
        ax[0].set_yticks(ticks=[0, 0.5, 1])
        ax[0].set_yticklabels(labels=[0, 0.5, 1], rotation=90, fontsize=8)
        ax[0].set_ylabel('SNP-index', fontsize=12, labelpad=15)
        ax[0].set_title('(a)', fontsize=17, rotation=0, x=-0.14, y=0.85)
        ax[0].axhline(y=0.5, color='black', linestyle='dashed', linewidth=0.75)
        ax[0].spines['top'].set_visible(False)
        ax[0].spines['right'].set_visible(False)

        # SNPidx2
        ax[1].scatter(x, y2, s=0.5, c=arg['--palette']['SNPidx2'], alpha=arg['--alpha'])
        ax[1].set(xlim=(0, max_x), ylim=(0, 1))
        ax[1].set_xticks(ticks=np.arange(0, max_x, 5e6))
        ax[1].tick_params(labelbottom=False)
        ax[1].set_yticks(ticks=[0, 0.5, 1])
        ax[1].set_yticklabels(labels=[0, 0.5, 1], fontsize=8, rotation=90)
        ax[1].set_ylabel('SNP-index', fontsize=12, labelpad=15)
        ax[1].set_title('(b)', fontsize=17, rotation=0, x=-0.14, y=0.85)
        ax[1].axhline(y=0.5, color='black', linestyle='dashed', linewidth=0.75)
        ax[1].spines['top'].set_visible(False)
        ax[1].spines['right'].set_visible(False)

        # D-SNPidx
        ax[2].scatter(x, y3, s=0.5, c=arg['--palette']['DELTA'], alpha=arg['--alpha'])
        ax[2].set(xlim=(0, max_x), ylim=(-1.5, 1.5))
        ax[2].set_yticks(ticks=[-1, 0, 1])
        ax[2].tick_params(labelbottom=True)
        ax[2].set_xticks(ticks=np.arange(0,max_x,5e6))
        ax[2].set_xticklabels(labels=np.arange(0, max_x, 5e6), fontsize=8)
        ax[2].set_xlabel(xlabel='Chromosomal position (bp)', fontsize=15)
        ax[2].xaxis.set_major_formatter(ScalarFormatter())
        ax[2].ticklabel_format(axis='x', style='scientific',scilimits=(6, 6), useMathText=True)
        ax[2].set_yticklabels(labels=[-1, 0, 1], rotation=90,fontsize=8)
        ax[2].set_ylabel('$\Delta$'+' (SNP-index)',fontsize=12, labelpad=15)
        ax[2].set_title('(c)', fontsize=17, rotation=0, x=-0.14, y=0.85)
        ax[2].axhline(y=0, color='black', linestyle='dashed', linewidth=0.75)
        ax[2].spines['top'].set_visible(False)
        ax[2].spines['right'].set_visible(False)

        if arg['--moving-avg'] != False:
            plot_avg_qtl_SNPidx(d, arg, ax)
            if arg['--ci95']:
                calc_ci(d, arg, ax[2])
        # Save file
        fig.subplots_adjust(hspace=0.1)
        rtch = rt.format(chrom[i])
        filename=rtch + typ
        filename=check_save(arg, filename)
        plt.savefig(arg['--outdir']+filename)
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
            write_caption(f,cap)

def calc_ci(d, arg, ax):
    z95=abs(st.norm.ppf(.025/arg['n_markers']))
    d['DP']=d['DPref_1']+d['DPalt_1']+d['DPref_2']+d['DPalt_2']
    d['DP1']=d['DPref_1']+d['DPalt_1']
    d['DP2']=d['DPref_2']+d['DPalt_2']
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

def plot_avg(d, arg, ax, field):
    if field == 'SNPidx2':
        d['DP']=d['DPref_2']+d['DPalt_2']
        c_=arg['--palette']['SNPidx2']
    elif field == 'SNPidx1':
        d['DP']=d['DPref_1']+d['DPalt_1']
        c_=arg['--palette']['SNPidx1']
    elif field == 'MAX_SNPidx2':
        if 'SNPidx2' not in arg['--fields']:
            d['DP']=d['DPref_1']+d['DPalt_1']
        else:
            d['DP']=d['DPref_2']+d['DPalt_2']
        c_=arg['--palette']['MAX_SNPidx2']
    elif field == 'log10PVALUE':
        d['DP']=d['DPref_1']+d['DPalt_1']+d['DPref_2']+d['DPalt_2']
        c_=arg['--palette']['log10PVALUE']
    elif field == 'ED100_4':
        d['DP']=d['DPref_1']+d['DPalt_1']+d['DPref_2']+d['DPalt_2']
        c_=arg['--palette']['log10PVALUE']
    elif field == 'DELTA':
        d['DP']=d['DPref_1']+d['DPalt_1']+d['DPref_2']+d['DPalt_2']
        c_=arg['--palette']['DELTA']
    elif field == 'G':
        d['DP']=d['DPref_1']+d['DPalt_1']+d['DPref_2']+d['DPalt_2']
        c_=arg['--palette']['mvg']
    d['prody']=d[field] * d['DP']
    d['avgy']=(d['prody'].rolling(arg['--moving-avg']).sum() / \
               d['DP'].rolling(arg['--moving-avg']).sum())
    d['prodx']=d['POS'] * d['DP']
    d['avgx']=(d['prodx'].rolling(arg['--moving-avg']).sum() / \
               d['DP'].rolling(arg['--moving-avg']).sum())
    ax.plot(d['avgx'], d['avgy'], c=c_, lw=1.25)#lw=0.75

def plot_avg_qtl_SNPidx(d, arg, ax):
    color = arg['--palette']['mvg']
    d['DP']=d['DPref_1']+d['DPalt_1']+d['DPref_2']+d['DPalt_2']
    d['DP1']=d['DPref_1']+d['DPalt_1']
    d['DP2']=d['DPref_2']+d['DPalt_2']
    d['productx']=d['POS'] * d['DP']
    d['mediamovilx']=(d['productx'].rolling(arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())
    #SNPidx1
    d['producty1']=d['SNPidx1'] * d['DP1']
    d['mediamovily1']=(d['producty1'].rolling(arg['--moving-avg']).sum() / d['DP1'].rolling(arg['--moving-avg']).sum())
    ax[0].plot(d['mediamovilx'], d['mediamovily1'], c=color, lw=1.5)
    #SNPidx2
    d['producty2']=d['SNPidx2'] * d['DP2']
    d['mediamovily2']=(d['producty2'].rolling(arg['--moving-avg']).sum() / d['DP2'].rolling(arg['--moving-avg']).sum())
    ax[1].plot(d['mediamovilx'], d['mediamovily2'], c=color, lw=1.5)
    # D-SNPidx
    d['productyD']=d['DELTA'] * d['DP']
    d['mediamovilD']=(d['productyD'].rolling(arg['--moving-avg']).sum() / d['DP'].rolling(arg['--moving-avg']).sum())
    ax[2].plot(d['mediamovilx'], d['mediamovilD'], c=color, lw=1.5)


def plot_boost(d, arg, ax, max_x):
    # Mediamovil Boost
    if 'SNPidx2' not in arg['--fields']:
        d['DP']=d['DPref_1']+d['DPalt_1']
    else:
        d['DP']=d['DPref_2']+d['DPalt_2']
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

def create_caption(arg, res_tit):
    f=open(arg['captions_dir'] + res_tit + '.txt', 'a')
    return f

def write_caption(f, text):
    for sentence in text:
        f.write(sentence+'\n')
    f.write('\n')
    f.close()
