from docopt import docopt
import sys
import os
import time
# sys.path.append(os.getcwd())
from pylibs.analysis_lib import *
from pylibs.graphics_lib import *


def mbs(argv):
    mbs_doc = """
    Usage:
       maptools.py mbs [options]
       maptools.py mbs

    Options:
      -h, --help                  Show this screen.
      -v, --version               Show the version
      -i, --input FILE            VCF input file. Can also come from a pipe.

    Input Options:
      -d, --data LIST             Pools genotype: dominant(D), recessive(R), parental dominant(Pd) and
                                  parental recessive(Pr) [default: D,R].
      -r, --ref-genotype STR      Which parental houses the reference, \"miss\" for missing genotype [default: D].
      -C, --max-depth INT         Maximum allele depth [default: 120].
      -c, --min-depth INT         Minimum allele depth [default: 20].
      -Q, --max-ratio INT         Maximum allele frequency in dominant pool [default: 100].
      -q, --min-ratio INT         Minimum allele frequency in dominant pool [default: 0].
      -M, --mutagen STR           Type of mutagen used to filter variants[default: EMS].

    Output Options:
      -o, --output FILE           Write output to file.
      -O, --output-type TYPE      \"txt\" tab separated, \"csv\" comma separated [default: txt]

    """

    arg = docopt(mbs_doc, argv=None, help=True,version=v_mbs)
    arg['pipe'] = sys.stdin.isatty()
    arg = check_args(mbs_doc, arg)
    arg['version'] = v_mbs
    output = arg['--output']
    fsal = False
    if output != None:
        fsal = open(arg['--output'], 'w')
    arg['fsal'] = fsal
    first = True
    if arg['--input']:
        inp = arg['inp']
    else:
        inp = sys.stdin
    choose_header(arg)
    write_argv(arg, argv)
    print(arg)
    for line in inp:
        if line.startswith('#'):
          read_header(arg,line)
        else:
            fields, pools = vcf_line_parser2(line, arg)
            if (fields, pools) != (0, 0):
                REF = fields[2]
                al_count = normalize(pools, REF, arg)
                if al_count != None:
                    flag = filter_mut(arg, al_count[:2])
                    if flag:
                      calcs = mbs_calc(al_count[2:], arg)
                      if calcs != None:
                          first = new_line(fsal, arg, first,
                                          fields[:2], al_count, calcs)


def qtl(argv):
    qtl_doc = """
  Usage:
     maptools.py qtl [options]
     maptools.py qtl

  Options:
    -h, --help                    Show this screen.
    -v, --version                 Show the version.
    -i, --input FILE              VCF input file. Can also come from a pipe.

  Input Options:
    -d, --data LIST               Pools genotype: dominant(D), recessive(R), parental dominant(Pd) and 
                                  parental recessive(Pr) [default: D,R].
    -r, --ref-genotype STR        Which parental houses the reference, \"miss\" for missing genotype [default: D].
    -C, --max-depth INT           Maximum read depth [default: 120].
    -c, --min-depth INT           Minimum read depth [default: 20].

  Output Options:
    -o, --output FILE             Write output file.
    -O, --output-type TYPE        \"txt\" tab separated, \"csv\" comma separated [default: txt]

  """
    arg = docopt(qtl_doc, argv=None, help=True, version=v_qtl)
    # print(arg)
    arg['pipe'] = sys.stdin.isatty()
    arg = check_args(qtl_doc, arg)
    arg['version'] = v_qtl
    output = arg['--output']
    fsal = False
    if output != None:
        fsal = open(arg['--output'], 'w')
    arg['fsal'] = fsal
    first = True
    if arg['--input']:
        inp = arg['inp']
    else:
        inp = sys.stdin
    choose_header(arg)
    write_argv(arg, argv)
    print(arg)
    for line in inp:
        if line.startswith('#'):
          read_header(arg,line)
        else:
            fields, pools = vcf_line_parser2(line, arg)
            if (fields, pools) != (0, 0):
                REF = fields[2]
                #print(fields[1])
                al_count = normalize(pools, REF, arg)
                if al_count != None:
                    calcs = qtl_calc(al_count[2:], arg)
                    if calcs != None:
                        first = new_line(fsal, arg, first,
                                         fields[:2], al_count, calcs)

def merge(argv):
  merge_doc="""
  Usage:
     maptools.py merge [options]
     maptools.py merge

  Options:
    -h, --help                    Show this screen.
    -v, --version                 Show the version.
    -i, --input FILE              Maptools MBS or QTL-Seq output.

  Input Options:                   
    -w, --window INT              Number of adjacent markers to merge the data [default: 20].
    -c, --chromosomes LIST        Chromosomes names selection (separeted by comma) [default: all].
  
  Output Options:
    -o, --output FILE             Write output to file.
    -O, --output-type TYPE        \"txt\" tab separated, \"csv\" comma separated [default: txt]
  """
  arg = docopt(merge_doc, argv=None, help=True, version=v_merge)
  arg = check_merge(arg, merge_doc)
  fsal = False
  if arg['--output'] != None:
    fsal = open(arg['--output'], 'w')
  arg['fsal'] = fsal
  arg['version'] = v_merge
  header_lines = read_header_merge(arg)
  write_argv(arg, argv)
  arg, df = load_dataframe(arg)
  #print(arg)
  for line in header_lines:
    write_line(line, fsal)
  grouped_by(df, arg)
  if fsal == True:
    fsal.close()

def mbs_plot(argv):
    mbsplot_doc = """
  Usage:
     maptools.py mbsplot <input_file> [options]
     maptools.py mbsplot

  Options:
    --chromosomes,-c=<opt>     List of chromosome names (separeted by comma) [default: all].
    --multi-chrom, -m          Multi-chromosome plots.
    --pvalue, -p               Generates p-value plots.
    --allele-freq-1, -D        Generates allele frequency plots for the dominant pool (AF1).
    --allele-freq-2, -R        Generates allele frequency plots for the recessive pool (AF2).
    --combine, -X              Combined plots with AF1 and AF2 lines. Use it together with --moving-avg.
    --max-allele-freq2, -M     Represents the maximum allele frequency in recessive pool. Recomended when data is not phased.
    --all, -a                  Generates all possible plot types.
    --moving-avg, -A=<int>     Number of adjacent markers for the calculation of moving average curves.
    --boost, -b=<int>          Number of adjacent markers for the calculation of moving average curves for boost, showed in recessive pool.
    --palette, -P=<int>        Select the colour pallete for your plots (1 for colour blindness) [default: 1].
    --alpha, -t=<float>        Select the dot transparency (0.0 to 1.0) [default: 0.4].
    --bonferroni               Show Bonferroni test line in p-value graphics.
    --fileformat=<type>        Output format: pdf, svg, jpg [default: pdf].
    --captions                 Generates figure captions.
    --outdir,-o=<dir>          Output directory [default: graphics].
  """
    arg = docopt(mbsplot_doc, argv=None, help=True,version=v_mbsplot)
    arg = test_plot(arg, mbsplot_doc)
    arg = read_header_plot(arg)
    arg, df = load_dataframe_plotting(arg)
    #print(arg)
    arg['version'] = v_mbsplot
    if arg['--pvalue'] == True:
        pval_mono_graph(df, arg)
        if arg['--multi-chrom'] == True:
          pval_multi_graph(df, arg)
          pval_multi_Vertical_graph(df, arg)
    if arg['--allele-freq-1'] == True and arg['--allele-freq-2'] == True:
        if arg['--combine'] and arg['--moving-avg'] != False:
            AF1_AF2_mono_graph(df, arg)
        if arg['--multi-chrom'] == True:
            AF_multi_Vertical_graph(df, arg, 'SNPidx1')
            AF_multi_Vertical_graph(df, arg, 'SNPidx2')
            if arg['--combine'] and arg['--moving-avg'] != False: 
              AF12_multi_Vertical_graph(df, arg)
    if arg['--allele-freq-1'] == True and arg['--allele-freq-2'] == True and arg['--pvalue'] == True:
        print('Allele Frequencies & P-value vertical single chromosome')
        AF1_AF2_pval_mono(df, arg)
    if arg['--allele-freq-1'] == True:
        print('Allele Frequency 1 single chromosome')
        AF_mono_graph(df, arg, 'SNPidx1')
        if arg['--multi-chrom'] == True:
           AF_multi_Vertical_graph(df, arg, 'SNPidx1')
    if arg['--allele-freq-2'] == True:
        print('Allele Frequency 2 single chromosome')
        AF_mono_graph(df, arg, 'SNPidx2')
        if arg['--multi-chrom'] == True:
           AF_multi_Vertical_graph(df, arg, 'SNPidx2')

    if arg['--max-allele-freq2'] == True:
        AF_mono_graph(df, arg, 'MAX_SNPidx2')
        if arg['--multi-chrom'] == True:
           AF_multi_Vertical_graph(df, arg, 'MAX_SNPidx2')


def qtl_plot(argv):
    qtlplot_doc = """QTL-Plot
  Usage:
     maptools.py qtlplot <input_file> [options]
     maptools.py qtlplot -v
     maptools.py qtlplot

  Options:
    -v --version               Show the version.

  Output options:
    --chromosomes,-c=<opt>     List of chromosome names (separeted by comma) [default: all].
    --multi-chrom, -m          Multi-chromosome plots.
    --pvalue, -p               Generates p-value plots.
    --delta, -d                Generates the SNPidx for high and low pools and Delta plot.
    --allele-freq-H, -H        Generates the alelle frequency plots for the pool of high phenotype.
    --allele-freq-L, -L        Generates the allele frecuency plots for the pool of low phenotype.
    --combine, -X              Combined plots with high and low lines. Use it together with --moving-avg.
    --euclidean-distance, -E   Generates Euclidean distance plots between individual SNPs and in groups of 100 adjacent markers.
    --g-statistic, -G          Generates G-statistic plots for individual SNPs.
    --qtl-seq, -Q              Generates muti-plots with ED, G, DELTA and p-value graphics for each chromosome.
    --all, -a                  Generates all possible plot types.
    --moving-avg, -A=<int>     Number of adjacent markers for the calculation of moving average curves.
    --palette, -P=<int>        Select the colour palette for your plots (1 for colour blidness) [default: 1].
    --alpha, -t=<float>        Select the dot transparency (0.0 to 1.0) [default: 0.4].
    --bonferroni               Show Bonferroni test line in p-value plots.
    --ci95                     Show the  95% cofindence interval in Delta plots. 
    --fileformat=<type>        Output format: pdf, svg, jpg [default: pdf].
    --captions, -C             Generates figure captions.
    --outdir,-o=<dir>          Output directory [default: graphics].

  """
    arg = docopt(qtlplot_doc, argv=None, help=True,version=v_qtlplot)
    arg = test_plot(arg, qtlplot_doc)
    arg = read_header_plot(arg)
    arg, df = load_dataframe_plotting(arg)
    arg['version'] = v_qtlplot
    #print(arg)
    if arg['--euclidean-distance']:
        df = get_ED100_4(df, arg, RANG)
        plot_ED(df,arg)
        if arg['--multi-chrom'] == True:
            ED_multi_Vertical_graph(df, arg)
    if arg['--g-statistic'] == True:
        plot_G(df, arg)
        if arg['--multi-chrom'] == True:
            G_multi_Vertical_graph(df, arg)
    if arg['--pvalue'] == True:
        pval_mono_graph(df, arg)
        if arg['--multi-chrom'] == True:    
            pval_multi_graph(df, arg)
            pval_multi_Vertical_graph(df, arg)
    if arg['--delta'] == True:
        snp_index_graph(df, arg)
        AF_mono_graph(df, arg, 'DELTA')
        if arg['--multi-chrom'] == True:
            AF_multi_Vertical_graph(df, arg, 'DELTA')
    
    if arg['--qtl-seq'] == True:
      df = get_ED100_4(df, arg, RANG)
      qtl_mixed_plot(df, arg)
    if arg['--allele-freq-H'] == True and arg['--allele-freq-L'] == True:
        if arg['--combine'] and arg['--moving-avg'] != False:
            AF1_AF2_mono_graph(df, arg)
        if arg['--multi-chrom'] == True:
            AF_multi_Vertical_graph(df, arg, 'SNPidx1')
            AF_multi_Vertical_graph(df, arg, 'SNPidx2')
            if arg['--combine'] and arg['--moving-avg'] != False:
                AF12_multi_Vertical_graph(df, arg)

    if arg['--allele-freq-H'] == True:
        AF_mono_graph(df, arg, 'SNPidx1')
    if arg['--allele-freq-L'] == True:
        AF_mono_graph(df, arg, 'SNPidx2')   
        
def annotate(argv):
  annotate_doc="""Annotate variants
Usage:
   maptools.py annotate [options]
   maptools.py annotate --version
   maptools.py annotate -h

Options:
   --help -h                     Show this screen.
   --version                     Show the version.
   --input, -i=<file>            VCF input file. If the input come from a pipe, don't use this option.
   --gff, -g=<file>              GFF3 input file.
   --reference, -f=<file>        FASTA reference file (The same used to align the reads).
Input Options:
   --data, -d=<opt>              Code of genotypes ordered according VCF input file [default: D,R].
   --ref, -r=<opt>               Which parental houses the reference [default: D].
   --no-ref                      Don't normalize data.
   --mutagen, -M=<opt>           Type of mutagen used to filter variants[defult: EMS].
   --region, -R=<region>         Region of the genome to explore (... -R chrName:Startpos-Endpos)
Output Options:
   --output, -o=<file>           Output file.
   --outdir, -O=<dir>            Output directory [default: results].           
  """
  arg = docopt(annotate_doc, argv=None, help=True, version=v_annotate)
  arg['pipe'] = sys.stdin.isatty()
  arg['version'] = v_annotate
  arg = test_arg_ann(annotate_doc, arg)
  df = create_df(arg)
  output = arg['--output']
  fsal = False
  arg['fsal'] = False
  if output != None:
    fsal = open(arg['--outdir']+arg['--output'], 'w')
    arg['fsal'] = fsal
  
  if arg['--input']:
    inp = arg['inp']
  else:
    inp = sys.stdin
  write_argv(arg, argv)
  for line in inp:
    if line.startswith('#'):
      read_header(arg,line)
    else:
      line = filter_region(line, arg)
      if line != None:
        fields, pools = vcf_line_parser2(line, arg)
        if (fields, pools) != (0,0):
          REF = fields[2]
          all_count = normalize_annotate(pools, REF, arg, fields)
          if all_count != None:
            flag = filter_mut(arg, all_count[:2])
            if flag:
              calcs = ann_calc(all_count[2:], arg)
              if calcs != None:
                df = new_df_line(df,arg,fields[:2], all_count, calcs)
  load_reference(df,arg)
  df = df.reset_index()
  start = time.perf_counter()
  load_gff(arg)
  write_annotate_header(arg)
  for idx, row in df.iterrows():
    check_mutation2(row,arg)
  finish = time.perf_counter()
  #print(f'Variant annotation finished in {round(finish-start, 2)} second(s)')
  
  ##
