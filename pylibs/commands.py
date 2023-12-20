from docopt import docopt
import sys
import os
import time
# sys.path.append(os.getcwd())
from pylibs.analysis_lib import *
from pylibs.graphics_lib import *
from pylibs.constants import *


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
  -d, --data LIST               comma separated list of sequenced samples
                                Available options:
                                  D   bulk of phenotypically dominant individuals from a segregating population
                                  R   bulk of phenotypically recessive individuals from a segregating population
                                  Pd  re-sequencing of the dominant parent genome
                                  Pr  re-sequencing of the recessive parent genome
                                  Wd  re-sequencing of the wild-type strain isogenic to a dominant mutant
                                  Wr  re-sequencing of the wild-type strain isogenic to a recessive mutant
                                    Note: Maptools requires that Wd, Wr, Pd and Pr are highly homozygous inbred lines
                              
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
    #print(arg)
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
            fields, pools, genotypes = vcf_line_parser3(line, arg) #nueva linea HC
            if (fields, pools, genotypes) != (0, 0, 0):
              #ref_allele = fields[2]
              #REPONER pools, genotypes, normalized = normalize2(pools, arg, genotypes) 
              normalized, sorted = normalize2(pools, arg, genotypes)
              ####print("Normalize2", pools, genotypes, normalized)
              #REPONERif (pools, genotypes, normalized) != (0,0,0):
              if normalized != False:
                ####print("Todo",pools,genotypes,normalized) #,al_normalized)
                if arg['--no-filter'] == False:
                   #flag = filter_mbs(arg, pools,genotypes,normalized,al_normalized)
                   flag = filter_mbs(arg, pools, genotypes, normalized, sorted)
                   ####print("Flag",flag)
                else:
                   flag = True
                if flag == True:
                   #print("Bien!!")
                   #print(al_normalized)
                   #print(al_normalized[2:])
                   calcs = mbs_calc(arg, pools, normalized, sorted) # antes era al_count
                   #print(calcs)
                   if calcs != None:
                      #al_normalized = 0 #pa que no de el error
                      #print(fields[:2], calcs)
                      alleles = [ allele for allele in pools['R'].keys()]
                      bases = [alleles[normalized[0]], alleles[normalized[1]]] 
                      first = new_line(fsal, arg, first, fields[:2] + bases, calcs)


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
                                    Note: Maptools requires that P samples are highly homozygous inbred lines

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
            fields, pools, genotypes = vcf_line_parser3(line, arg) #nueva linea HC
            if (fields, pools, genotypes) != (0, 0, 0):
              #ref_allele = fields[2]
              #REPONER pools, genotypes, normalized = normalize2(pools, arg, genotypes) 
              normalized, sorted = normalize2(pools, arg, genotypes)
              ####print("Normalize2", pools, genotypes, normalized)
              #REPONERif (pools, genotypes, normalized) != (0,0,0):
              if normalized != False:
                ####print("Todo",pools,genotypes,normalized) #,al_normalized)
                if arg['--no-filter'] == False:
                   #flag = filter_mbs(arg, pools,genotypes,normalized,al_normalized)
                   flag = filter_qtl(arg, pools, genotypes, normalized, sorted)
                   ####print("Flag",flag)
                else:
                   flag = True
                if flag == True:
                   #print("Bien!!")
                   #print(al_normalized)
                   #print(al_normalized[2:])
                   calcs = qtl_calc(arg, pools, normalized, sorted) # antes era al_count
                   #print(calcs)
                   if calcs != None:
                      #al_normalized = 0 #pa que no de el error
                      #print(fields[:2], calcs)
                      alleles = [ allele for allele in pools['R'].keys()]
                      bases = [alleles[normalized[0]], alleles[normalized[1]]] 
                      first = new_line(fsal, arg, first, fields[:2] + bases, calcs)

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
  arg, df = load_dataframe_merge(arg)
  #print(arg)
  for line in header_lines:
    write_line(line, fsal)
  grouped_by(df, arg)
  if fsal == True:
    fsal.close()

def plot(argv):
  plot_doc = """
Usage:
  maptools.py plot [options]
  maptools.py plot

Options:
  -h, --help                    show this help message and exit
  -v, --version                 print version information and exit
  -i, --input FILE              input file produced by the mbs or qtl commands
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
  --ci95                        add 95% confidence interval to delta plots. Requires -A or -W
  --palette STR                 select a colour palette [default: standard]
                                Available options: \"standard\", \"color_blind\" or \"custom\"
  -D, --dot-size FLOAT          Size of the points plotted [default: 1.0]
  -L, --line-width FLOAT        Width of trend lines [default: 2.0]
  --DPI INT                     Dots per inch [default: 600]

Plot types:
  -p, --pvalue                  generate p-value plots
  -1, --allele-freq-1           generate allele frequency plots for the D/high bulk (AF1)
  -2, --allele-freq-2           generate allele frequency plots for the R/low bulk (AF2)
  -M, --max-allele-freq2        plot the frequency of the most abundant allele in R bulk. Use -M when the
                                alleles cannot be assigned to the parents with certainty. Only for mbs analysis
  -d, --delta                   generate Delta plots (AF1 - AF2)
  -X, --combine                 overlay the moving average of AF2 on the AF1 plot, and vice versa. Requires
                                --moving-avg or --distance-avg.
  -E, --euclidean-distance      plot Euclidean distances, both individually (EDm) and as ED100^4 values
  -G, --g-statistic             generate G-statistic plots for individual markers
  -Q, --comb-statistics         generate multipanel figure with allele frequencies, DELTA, ED, G and p-value plots
                                for individual chromosomes
  -m, --multi-chrom             generate multi-chromosome plots for each statistic and the chromosomes 
                                especified with -c.
  -a, --all                     generate all possible plots
  
Output options:
  -o, --outdir DIR              output plot files to DIR [default: graphics]
  -O, --output-type TYPE        available types: pdf, svg, jpg [default: pdf]
  """
  arg = docopt(plot_doc, argv=None, help=True,version=v_plot)
  #print(arg)
  arg = test_plot(arg, plot_doc)
  arg = read_header_plot(arg)
  arg, df = load_dataframe_plotting(arg)

  #cc = arg['--contigs']
  #del arg['--contigs']
  #print(arg)
  #arg['--contigs'] = cc
  arg['version'] = v_plot
  #print(arg)
  #Allele Frequency plots
  if arg['--allele-freq-1'] == True:
    AF_mono_graph(df, arg, 'SNPidx1')
    if arg['--multi-chrom'] != False:
      AF_multi_Vertical_graph(df, arg, 'SNPidx1')
      AF_manhattan_plot(df, arg, 'SNPidx1')

  if arg['--allele-freq-2'] == True:
    AF_mono_graph(df, arg, 'SNPidx2')
    if arg['--multi-chrom'] != False:
      AF_multi_Vertical_graph(df, arg, 'SNPidx2')
      AF_manhattan_plot(df, arg, 'SNPidx2') 
  if arg['--combine'] == True:
    AF1_AF2_mono_graph(df, arg)
    if arg['--multi-chrom'] != False:
      AF12_multi_Vertical_graph(df, arg)
      AFCombinedManhattanPlot(df, arg)
  
  if arg['--max-allele-freq2'] == True and arg['type'] == 'mbs':
    AF_mono_graph(df, arg, 'MAX_SNPidx2')
    if arg['--multi-chrom'] != False:
      AF_multi_Vertical_graph(df, arg, 'MAX_SNPidx2')
      AF_manhattan_plot(df, arg, 'MAX_SNPidx2')
  #Delta plots
  if arg['--delta'] == True:
    AF_mono_graph(df, arg, 'DELTA')
    if arg['--multi-chrom'] == True:
      AF_multi_Vertical_graph(df, arg, 'DELTA')
      AF_manhattan_plot(df, arg, 'DELTA')
  
  #P-value plots
  if arg['--pvalue'] == True:
      pval_mono_graph(df, arg)
      if arg['--multi-chrom'] == True:
        #pval_multi_graph(df, arg)
        pval_multi_Vertical_graph(df, arg)
        pval_manhattan_plot(df, arg)
  ###REVISAR
  #Euclidean Distance plots
  if arg['--euclidean-distance'] == True:
    EDmonoPlot(df,arg)
    if arg['--multi-chrom'] == True:
      ED_multi_Vertical_graph(df, arg)
      EDmanhattanPlot(df, arg)
  #G-statistic
  if arg['--g-statistic'] == True:
    GmonoPlot(df, arg)
    if arg['--multi-chrom'] == True:
      G_multi_Vertical_graph(df, arg)
      GmanhattanPlot(df,arg)
  
  if arg['--comb-statistics'] == True:
    combinedPlot(df, arg)
       

def annotate(argv):
  annotate_doc="""
Usage:
  maptools.py annotate [options]
  maptools.py annotate

Options:
  -h, --help                    show this help message and exit
  -v, --version                 print version information and exit
  -i, --input FILE              input file, mbs command output
  -g, --gff FILE                genome annotation file in gff3 format
  -f, --fasta-reference FILE    genome sequence file in fasta format

Input Options:
  -R, --region REGION           region of the genome to explore (-R chrName:Startpos-Endpos)
  -t, --transl-table INT        translation table [default: 1]

Output Options:
  -o, --output FILE             write output to FILE [standard output]
  -O, --output-type TYPE        txt: tab separated, csv: comma separated [default: txt]                     
  """
  arg = docopt(annotate_doc, argv=None, help=True, version=v_annotate)
  arg['pipe'] = sys.stdin.isatty()
  arg['version'] = v_annotate
  arg = check_args(annotate_doc, arg)
  #df = create_df(arg)
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
    arg['spacerin'] = '\t'
  write_argv(arg, argv)
  #print(arg)

  for line in inp:
    if line.startswith('#'):
      if line.startswith('##maptools_mbsCommand'):
        arg['mbs'] = True
        write_line(line,fsal)
        df = getAnalysisArgs(arg,line)
      if line.startswith('##maptools_qtlCommand'):
        arg['qtl'] = True
        write_line(line,fsal)
        df = getAnalysisArgs(arg,line) 
      read_header(arg,line)
    else:
      line = filter_region(line,arg)
      if line != None:
        fields = analysisParser(line, arg)
        df = newDFline2(df, arg, fields)
  if df.empty:
     print('Warning: There is no variants to analyse. Please reduce filtering.', file=sys.stderr)
     sys.exit()
  load_reference(df,arg)
  df = df.reset_index()
  load_gff(arg)
  for idx, row in df.iterrows():
    if len(row['DOM']) > 1 or len(row['REC']) > 1:
       is_indel(row, arg)
       continue
    arg['variant'] = 'substitution'
    check_mutation2(row,arg)


def annotatevcf(argv):
  annotate_doc="""
Usage:
  maptools.py annotatevcf [options]
  maptools.py annotatevcf

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
                                    Note: Maptools requires that Wd, Wr, Pd and Pr are highly homozygous inbred lines

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
      line = filter_region(line,arg)
      if line != None:
        #fields, pools, genotype = vcf_line_parser2(line, arg)
        fields, pools, genotypes = vcf_line_parser3(line, arg) #nueva linea HC
        if (fields, pools, genotypes) != (0, 0, 0):
           normalized, sorted = normalize2(pools, arg, genotypes)
           if normalized != False:
             if arg['--no-filter'] == False:
               flag = filter_mbs(arg, pools, genotypes, normalized, sorted)
             else:
               flag = True
             if flag == True:
               calcs = mbs_calc(arg, pools, genotypes, normalized, sorted) # antes era al_count
               if calcs != None:
                 alleles = [ allele for allele in pools['R'].keys()]
                 bases = [alleles[normalized[0]], alleles[normalized[1]]] 
#                first = new_line(fsal, arg, first, fields[:2] + bases, calcs)
                 #df = new_df_line(df,arg,fields[:2], al_count,calcs)
                 if normalized == {0:0, 1:1}:
                       reorder = [0]
                 else:
                       reorder = [1]
                 df = new_df_line(df,arg,fields[:2], bases, calcs,reorder)
  if df.empty:
     print('Warning: There is no variants to analyse. Please reduce filtering.', file=sys.stderr)
     sys.exit()
  load_reference(df,arg)
  df = df.reset_index()
  load_gff(arg)
  write_annotate_header(arg)
  for idx, row in df.iterrows():
    if len(row['DOM']) > 1 or len(row['REC']) > 1:
       is_indel(row, arg)
       continue
    arg['variant'] = 'substitution'
    check_mutation2(row,arg)
  
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
