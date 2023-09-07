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
  -h, --help                    show this help message and exit
  -v, --version                 print version information and exit
  -i, --input FILE              input file in uncompressed vcf format

Input Options:
  -d, --data LIST               comma separated list of sequenced samples [default: D,R]
                                Available options:
                                  D   bulk of phenotypically dominant individuals from a segregating population
                                  R   bulk of phenotypically recessive individuals from a segregating population
                                  Pd  re-sequencing of the dominant parent genome
                                  Pr  re-sequencing of the recessive parent genome
                                  Wd  re-sequencing of the wild-type strain isogenic to a dominant mutant
                                  Wr  re-sequencing of the wild-type strain isogenic to a recessive mutant
                              
  -r, --ref-genotype STR        indicate if the reference genome sequence corresponds to one of the parents
                                of the mapping population [default: miss]
                                Available options:
                                  D     reference genome corresponds to the phenotypically dominant parent
                                  R     reference genome corresponds to the phenotypically recessive parent
                                  miss  reference genome does not match the sequence of either parent

  -m, --mutant-pool STR         indicate the bulk established from phenotypically individuals from the
                                segregating population [default: R]
                                Available options:
                                  D   for a dominant mutation
                                  R   for a recessive mutation
Output Options:
  -o, --output FILE             write output to FILE [standard output]
  -O, --output-type TYPE        txt: tab separated, csv: comma separated [default: txt]
    
Filter Options:
  -C, --max-depth INT           maximum read depth in the D and R bulks for a position to be considered [default: inf]
  -c, --min-depth INT           minimum read depth in the D and R bulks for a position to be considered [default: 0]
  -Q, --max-ratio INT           maximum allele frequency in the D bulk [default: 100]
  -q, --min-ratio INT           minimum allele frequency in the D bulk [default: 0]
                                Note: you can set different -q and -Q values to enforce that markers are heterozygous
                                in the D bulk, by only considering those whose allele frequency in the interval [-q, -Q].
                                Requires --het-filter
                  
  --EMS                         ignore SNPs not caused by EMS (keeps G/C-to-A/T transition mutations)
  -I, --skip-indels             ignore indels
  --parental-filter             ignore variants shared with parental strain (requires one of: Pr, Pd,
                                Wr or Wd, specified with -d)
  --het-filter                  selects markers that are heterozygous in the D pool, either because the allele frequencies
                                are in the [-q, -Q] interval or because their GT is \"0/1\" in the vcf input
  --no-filter                   disable all filters              
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
    #print(arg)
    for line in inp:
        if line.startswith('#'):
          read_header(arg,line)
        else:
            fields, pools, genotype = vcf_line_parser2(line, arg)
            if (fields, pools, genotype) != (0, 0, 0):
              DOM = fields[2]
              arg['poss'] = fields[1]
              al_count,p_al_count,genotype = normalize(pools, DOM, arg, genotype)
              if (al_count,p_al_count,genotype) != (0,0,0):
                if arg['--no-filter'] == False:
                  flag = filter_mbs(arg,al_count,p_al_count, genotype)
                else:
                   flag = True
                if flag == True:
                  calcs = mbs_calc(al_count[2:], arg)
                  #print(calcs)
                  if calcs != None:
                    first = new_line(fsal, arg, first, fields[:2], al_count, calcs)


def qtl(argv):
    qtl_doc = """
Usage:
  maptools.py qtl [options]
  maptools.py qtl

Options:
  -h, --help                    show this help message and exit
  -v, --version                 print version information and exit
  -i, --input FILE              input file in uncompressed vcf format

Input Options:
  -d, --data LIST               comma separated list of sequenced samples [default: H,L]
                                Available options:
                                  H   bulk of individuals with the highest values for a quantitative trait
                                  L   bulk of individuals with the lowest values for a quantitative trait
                                  P   re-sequencing of a parent of the mapping population

  -r, --ref-genotype STR        indicate if the reference genome sequence corresponds to one of the parents
                                of the mapping population [default: miss]
                                Available options:
                                  P     reference genome corresponds to one of the parents of the cross
                                  miss  reference genome does not match the sequence of either parent
  
Filter Options:
  -C, --max-depth INT           maximum read depth in the H and L bulks for a position to be considered [default: inf]
  -c, --min-depth INT           minimum read depth in the H and L bulks for a position to be considered [default: 0]
  -I, --skip-indels             ignore indels
  --no-filter                   disable all filters

Output Options:
  -o, --output FILE             write output to FILE [standard output]
  -O, --output-type TYPE        txt: tab separated, csv: comma separated [default: txt]
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
    #print(arg)
    for line in inp:
        if line.startswith('#'):
          read_header(arg,line)
        else:
            fields, pools, genotype = vcf_line_parser2(line, arg)
            if (fields, pools, genotype) != (0, 0, 0):
                DOM = fields[2]
                al_count,p_al_count,genotype = normalize(pools, DOM, arg, genotype)
                if (al_count,p_al_count,genotype) != (0,0,0):
                    if arg['--no-filter'] == False:
                       flag = filter_qtl(arg, al_count, genotype)
                    else:
                       flag = True
                    if flag:
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
  -h, --help                    show this help message and exit
  -v, --version                 print version information and exit
  -i, --input FILE              input file produced by the mbs or qtl commands

Input Options:                   
  -w, --window INT              number of markers per bin
  -D, --distance INT            chromosome distance to cluster markers within a INT(bp)-range
  -c, --chromosomes LIST        comma-separated list of chromosome names [default: all]
  
Output Options:
  -o, --output FILE             write output to FILE [standard output]
  -O, --output-type TYPE        txt: tab separated, csv: comma separated [default: txt]
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
  maptools.py mbsplot [options]
  maptools.py mbsplot

Options:
  -h, --help                    show this help message and exit
  -v, --version                 print version information and exit
  -i, --input FILE              input file produced by the mbs command
  -c, --chromosomes LIST        comma-separated list of chromosome names [default: all]
  -C, --captions                automatically generate figure legends

Plot options:
  -A, --moving-avg INT          add moving averages to plots, calculated using INT adjacent markers
  -W, --distance-avg INT        add weighted averages to plots, calculated for all markers located
                                within an INT bp interval
  -b, --boost INT               add boost to allele frequency plots, calculated as the average boost
                                values of INT adjacent markers
  -B, --distance-boost INT      add boost to allele frequency plots, calculated for all markers located
                                within an INT bp interval
  -t, --alpha FLOAT             marker transparency in plots (0 - 1) [default: 0.4]
  --bonferroni                  add Bonferroni threshold to p-value plots. Requires -p or -a
  --palette STR                 select a colour palette [default: standard]
                                Available options: \"standard\", \"color_blind\" or \"custom\" 

Plot types:
  -p, --pvalue                  generate p-value plots
  -D, --allele-freq-1           generate allele frequency plots for the D bulk (AF1)
  -R, --allele-freq-2           generate allele frequency plots for the R bulk (AF2)
  -M, --max-allele-freq2        plot the frequency of the most abundant allele in R bulk. Use -M when the
                                alleles cannot be assigned to the parents with certainty
  -X, --combine                 overlay the moving average of AF2 on the AF1 plot, and vice versa. Requires
                                --moving-avg or --distance-avg.
  -m, --multi-chrom             generate multi-chromosome plots for each statistic and the chromosomes 
                                especified with -c. Manhattan plots require -p or -a
  -a, --all                     generate all possible plots
  
Output options:
  -o, --outdir DIR              output plot files to DIR [default: graphics]
  -O, --output-type TYPE        available types: pdf, svg, jpg [default: pdf]
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
          #pval_multi_graph(df, arg)
          pval_multi_Vertical_graph(df, arg)
          pval_manhattan_plot(df, arg)
    if arg['--allele-freq-1'] == True and arg['--allele-freq-2'] == True:
        if arg['--combine'] and (arg['--moving-avg'] != False or arg['--distance-avg'] != False):
            AF1_AF2_mono_graph(df, arg)
            if arg['--multi-chrom'] != False: 
              AF12_multi_Vertical_graph(df, arg)
    if arg['--allele-freq-1'] == True and arg['--allele-freq-2'] == True and arg['--pvalue'] == True:
        #print('Allele Frequencies & P-value vertical single chromosome')
        AF1_AF2_pval_mono(df, arg)
    if arg['--allele-freq-1'] == True:
        #print('Allele Frequency 1 single chromosome')
        AF_mono_graph(df, arg, 'SNPidx1')
        if arg['--multi-chrom'] == True:
           AF_multi_Vertical_graph(df, arg, 'SNPidx1')
    if arg['--allele-freq-2'] == True:
        #print('Allele Frequency 2 single chromosome')
        AF_mono_graph(df, arg, 'SNPidx2')
        if arg['--multi-chrom'] == True:
           AF_multi_Vertical_graph(df, arg, 'SNPidx2')

    if arg['--max-allele-freq2'] == True:
        AF_mono_graph(df, arg, 'MAX_SNPidx2')
        if arg['--multi-chrom'] == True:
           AF_multi_Vertical_graph(df, arg, 'MAX_SNPidx2')


def qtl_plot(argv):
    qtlplot_doc = """
Usage:
  maptools.py qtlplot [options]
  maptools.py qtlplot

Options:
  -h, --help                    show this help message and exit
  -v, --version                 print version information and exit
  -i, --input FILE              input file produced by the qtl command
  -c, --chromosomes LIST        comma-separated list of chromosome names [default: all]
  -C, --captions                automatically generate figure legends

Plot options:
  -A, --moving-avg INT          add moving averages to plots, calculated using INT adjacent markers
  -W, --distance-avg INT        add weighted averages to plots, calculated for all markers located
                                within an INT bp interval
  -t, --alpha FLOAT             marker transparency in plots (0 - 1) [default: 0.4]
  --bonferroni                  add Bonferroni threshold to p-value plots. Requires -p or -a
  --ci95                        add 95% confidence interval to delta plots. Requires -A or -W
  --palette STR                 select a colour palette [default: standard]
                                Available options: \"standard\", \"color_blind\" or \"custom\" 

Plot types:
  -p, --pvalue                  generate p-value plots
  -H, --allele-freq-H           generate allele frequency plots for the H bulk (AF1)
  -L, --allele-freq-L           generate allele frequency plots for the L bulk (AF2)
  -d, --delta                   generate Delta plots (AF1 - AF2)
  -X, --combine                 overlay the moving average of AF2 on the AF1 plot, and vice versa. Requires
                                --moving-avg
  -E, --euclidean-distance      plot Euclidean distances, both individually (EDm) and as ED100^4 values
  -G, --g-statistic             generate G-statistic plots for individual markers
  -Q, --qtl-seq                 generate multipanel figure with ED, G, DELTA and p-value plots
                                for individual chromosomes
  -m, --multi-chrom             generate multi-chromosome plots for each statistic and the chromosomes 
                                especified with -c. Manhattan plots require -p or -a 
  -a, --all                     generate all possible plots
  
Output options:
  -o, --outdir DIR              output plot files to DIR [default: graphics]
  -O, --output-type TYPE        available types: pdf, svg, jpg [default: pdf]
  """
    arg = docopt(qtlplot_doc, argv=None, help=True,version=v_qtlplot)
    #print(arg)
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
            #pval_multi_graph(df, arg)
            pval_multi_Vertical_graph(df, arg)
            pval_manhattan_plot(df, arg)
    if arg['--delta'] == True and arg['--ref-genotype'] == True:
        snp_index_graph(df, arg)
    
    if arg['--delta'] == True:
        AF_mono_graph(df, arg, 'DELTA')
        if arg['--multi-chrom'] == True:
            AF_multi_Vertical_graph(df, arg, 'DELTA')
    
    if arg['--qtl-seq'] == True:
      df = get_ED100_4(df, arg, RANG)
      qtl_mixed_plot(df, arg)
    if arg['--allele-freq-H'] == True and arg['--allele-freq-L'] == True:
        if arg['--combine'] and (arg['--moving-avg'] != False or arg['--distance-avg'] != False):
            AF1_AF2_mono_graph(df, arg)
        if arg['--multi-chrom'] == True:
            AF_multi_Vertical_graph(df, arg, 'SNPidx1')
            AF_multi_Vertical_graph(df, arg, 'SNPidx2')
            if arg['--combine'] and (arg['--moving-avg'] != False or arg['--distance-avg'] != False):
                AF12_multi_Vertical_graph(df, arg)

    if arg['--allele-freq-H'] == True:
        AF_mono_graph(df, arg, 'SNPidx1')
    if arg['--allele-freq-L'] == True:
        AF_mono_graph(df, arg, 'SNPidx2')   
        
def annotate(argv):
  annotate_doc="""
Usage:
  maptools.py annotate [options]
  maptools.py annotate

Options:
  -h, --help                    show this help message and exit
  -v, --version                 print version information and exit
  -i, --input FILE              input file in uncompressed vcf format
  -g, --gff FILE                genome annotation file in gff3 format
  -f, --fasta-reference FILE    genome sequence file in fasta format

Input Options:
  -d, --data LIST               comma separated list of sequenced samples [default: D,R]
                                Available options:
                                  D   bulk of phenotypically dominant individuals from a segregating population
                                  R   bulk of phenotypically recessive individuals from a segregating population
                                  Pd  re-sequencing of the dominant parent genome
                                  Pr  re-sequencing of the recessive parent genome
                                  Wd  re-sequencing of the wild-type strain isogenic to a dominant mutant
                                  Wr  re-sequencing of the wild-type strain isogenic to a recessive mutant

  -r, --ref-genotype STR        indicate if the reference genome sequence corresponds to one of the parents
                                of the mapping population [default: miss]
                                Available options:
                                  D     reference genome corresponds to the phenotypically dominant parent
                                  R     reference genome corresponds to the phenotypically recessive parent
                                  miss  reference genome does not match the sequence of either parent

  -m, --mutant-pool STR         indicate the bulk established from phenotypically individuals from the
                                segregating population [default: R]
                                Available options:
                                  D   for a dominant mutation
                                  R   for a recessive mutation

  -R, --region REGION           region of the genome to explore (-R chrName:Startpos-Endpos)
  -t, --transl-table INT        translation table [default: 1]

Output Options:
  -o, --output FILE             write output to FILE [standard output]
  -O, --output-type TYPE        txt: tab separated, csv: comma separated [default: txt]

Filter Options:
  -C, --max-depth INT           maximum read depth in the D and R bulks for a position to be considered [default: inf]
  -c, --min-depth INT           minimum read depth in the D and R bulks for a position to be considered [default: 0]
  -Q, --max-ratio INT           maximum allele frequency in the D bulk [default: 100]
  -q, --min-ratio INT           minimum allele frequency in the D bulk [default: 0]
                                Note: you can set different -q and -Q values to enforce that markers are heterozygous
                                in the D bulk, by only considering those whose allele frequency in the interval [-q, -Q].
                                Requires --het-filter
                  
  --EMS                         ignore SNPs not caused by EMS (keeps G/C-to-A/T transition mutations)
  -I, --skip-indels             ignore indels
  --parental-filter             ignore variants shared with parental strain (requires one of: Pr, Pd,
                                Wr or Wd, specified with -d)
  --het-filter                  selects markers that are heterozygous in the D pool, either because the allele frequencies
                                are in the [-q, -Q] interval or because their GT is \"0/1\" in the vcf input
  --no-filter                   disable all filters                         
  """
  arg = docopt(annotate_doc, argv=None, help=True, version=v_annotate)
  arg['pipe'] = sys.stdin.isatty()
  arg['version'] = v_annotate
  arg['mbs'] = True
  arg = check_args(annotate_doc, arg)
  df = create_df(arg)
  output = arg['--output']
  fsal = False
  arg['fsal'] = False
  if output != None:
    fsal = open(arg['--output'], 'w')
  arg['fsal'] = fsal
  
  if arg['--input']:
    inp = arg['inp']
  else:
    inp = sys.stdin
  write_argv(arg, argv)
  #print(arg)
  for line in inp:
    if line.startswith('#'):
      read_header(arg,line)
    else:
      line = filter_region(line, arg)
      if line != None:
        fields, pools, genotype = vcf_line_parser2(line, arg)
        if (fields, pools, genotype) != (0,0,0):
          DOM = fields[2]
          al_count,p_al_count,genotype = normalize(pools, DOM, arg, genotype)
          if (al_count,p_al_count,genotype) != (0,0,0):
            if arg['--no-filter'] == False:
              flag = filter_mbs(arg, al_count, p_al_count, genotype)
            else:
              flag = True
            if flag:
              calcs = mbs_calc(al_count[2:], arg)
              if calcs != None:
                df = new_df_line(df,arg,fields[:2], al_count, calcs)
  if df.empty:
     print('Warning: There is no variants to analyse. Please reduce filtering.', file=sys.stderr)
     sys.exit()
  load_reference(df,arg)
  #print(str(arg['ref'].seq[81033575-1:81033653+5]))
  df = df.reset_index()
  start = time.perf_counter()
  load_gff(arg)
  write_annotate_header(arg)
  for idx, row in df.iterrows():
    if len(row['DOM']) > 1 or len(row['REC']) > 1:

       is_indel(row, arg)
       continue
    arg['variant'] = 'substitution'
    check_mutation2(row,arg)
  finish = time.perf_counter()
  #print(f'Variant annotation finished in {round(finish-start, 2)} second(s)')
  
  ##
def citation(argv):
  cite_doc="""
Program: MAPtools
Version: 1.00
Contact: hcandela@umh.es

If you find MAPtools useful, please include the following citation in your scientific publications:

Martínez-Guardiola, C., Parreño, R. and Candela, H. (2023). MAPtools: Command-Line Tools for
Mapping-by-Sequencing and QTL-Seq Analysis and Visualization. Submitted.                     
  """
  print(cite_doc, file=sys.stderr)