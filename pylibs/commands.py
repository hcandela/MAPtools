from docopt import docopt
import sys
import os
import time
# sys.path.append(os.getcwd())
from pylibs.analysis_lib import *
from pylibs.graphics_lib import *


def mbs(argv):
    mbs_doc = """Mapping by Sequencing analysis
    Usage:
       maptools.py mbs [options]
       maptools.py mbs --version
       maptools.py mbs -h

    Options:
      --help -h                   Show this screen.
      --version                   Show the version
      --input, -i=<file>          VCF input file. Can algo come from a pipe.
    Input Options:
      --data,-d=<opt>             Code of Genotypes [default: D,R].
      --ref, -r=<opt2>            Which parental houses the reference [default: D].
      --no-ref                    Don't normalize data
      --max-depth=<int>           Filt out reads with depth above maximum, in dominant pool [default: 120]
      --min-depth=<int>           Filt out reads with depth below minimum, in dominant pool [default: 20]
      --max-ratio=<int>           Filt out reads with allele frequency ratio above maximum, in dominant pool [default: 85]
      --min-ratio=<int>           Filt out reads with allele frequency ratio below minimum, in dominant pool [default: 15]
    Output Options:
      --output, -o=<file>         Output file
      --outdir, -O=<dir>          Output directory [default: results]
      --fileformat=<ext>          Formats availables: txt, csv (,) [default: txt]
    """

    arg = docopt(mbs_doc, argv=None, help=True,
                 version='Mapping by Sequencing analysis version: 0.1')
    arg['pipe'] = sys.stdin.isatty()
    arg = test_args(mbs_doc, arg)

    output = arg['--output']
    fsal = False
    if output != None:
        fsal = open(arg['--outdir']+arg['--output'], 'w')
    first = True
    if arg['--input']:
        inp = open(arg['--input'], 'r')
    else:
        inp = sys.stdin
    for line in inp:
        if not line.startswith('#'):
            fields, pools = vcf_line_parser2(line, arg)
            if (fields, pools) != (0, 0):
                REF = fields[2]
                al_count = normalize(pools, REF, arg)
                if al_count != None:
                    calcs = mbs_calc(al_count[2:], arg)
                    if calcs != None:
                        first = new_line(fsal, arg, first,
                                         fields[:2], al_count, calcs)


def qtl(argv):
    qtl_doc = """QTL-Seq analysis
  Usage:
     maptools.py qtl [options]
     maptools.py qtl --version
     maptools.py qtl -h

  Options:
    --help -h                   Show this screen.
    --version                   Show the version.
    --input, -i=<file>          VCF input file. Can algo come from a pipe.
  Input Options:
    --data,-d=<opt>             Code of Genotypes [default: D,R].
    --ref, -r=<opt2>            Linear references used in mapping step [default: D].
    --no-ref                    Don't normalize data
    --max-depth=<int>           Filt out positions with total (High + Low) depth above maximum [default: 120]
    --min-depth=<int>           Filt out positions with total (High + Low) depth below minimum [default: 20]
  Output Options:
    --output, -o=<file>         Output file
    --outdir, -O=<dir>          Output directory [default: results]
    --fileformat=<ext>          Formats availables: txt (tab), csv (,) [default: txt]

  """
    arg = docopt(qtl_doc, argv=None, help=True, version='0.1')
    # print(arg)
    arg['pipe'] = sys.stdin.isatty()
    arg = test_args(qtl_doc, arg)
    output = arg['--output']
    fsal = False
    if output != None:
        fsal = open(arg['--outdir']+arg['--output'], 'w')
    first = True
    if arg['--input']:
        inp = open(arg['--input'], 'r')
    else:
        inp = sys.stdin
    for line in inp:
        if not line.startswith('#'):
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
  merge_doc="""Merge
  Usage:
     maptools.py merge <input_file> [options]
     maptools.py merge -v
     maptools.py merge -h

  Options:
     -h --help                     Show this screen.
     -v --version                  Show the version.
     --window, -w=<int>            Number of adjacent markers to merge the data [default: 20].
     --chromosomes, -c=<list_int>  List of chromosomes separated by comma, from 1 to nº chr [default: all].
     --outdir, -O=<dir>            Create a new directory if it doesn't exist. If not provide uses the working directory [default: merge]
     --output, -o=<file>           Output file.
     --fileformat=<ext>            Formats availables: txt (tab), csv (,) [default: txt]
  """
  arg = docopt(merge_doc, argv=None, help=True, version='Merge mbs or qtl analysis version=0.1')
  df, arg = test_merge(arg)
  fsal = False
  if arg['--output'] != None:
    fsal = open(arg['--outdir']+arg['--output'], 'w')
  arg['fsal'] = fsal
  grouped_by(df, arg)
  if fsal == True:
    fsal.close()

def mbs_plot(argv):
    mbsplot_doc = """MBS
  Usage:
     maptools.py mbsplot <input_file> [options]
     maptools.py mbsplot -v
     maptools.py mbsplot -h

  Options:
    -h --help                  Show this screen.
    -v --version               Show the version.

  Output options:
    --chromosomes,-c=<opt>     List of chromosomes separeted by comma, from 1 to nºchr [default: all].
    --multi-chrom, -m          Multigraph representation.
    --pvalue, -p               Generates p-value graphics.
    --allele-freq-1, -D        Generates the allele frequency graphics for the recessive pool (AF1).
    --allele-freq-2, -R        Generates allele frequency graphics for the dominant pool (AF2).
    --combine, -C              Combined multi-graphics with AF1 and AF2 lines. Use it together with --moving-avg.
    --max-allele-freq2, -M     Represents the maximum allele frequency in recessive pool. Recomended when data is not phased.
    --all, -a                  Generates all possible graphics.
    --moving-avg, -A=<int>     Selects the number of adjacent markers for the construction of moving average curves for allele frequncy and p-value.
    --boost, -b=<int>          Selects the number of adjacent markers for the construction of moving average curves for boost, showed in recessive pool.
    --palette, -P=<int>        Select the colour pallete for your grpahics (1 for colour blindness) [default: 1].
    --alpha, -t=<float>        Select the dots transparency on graphics (0.0 to 1.0) [default: 0.4].
    --bonferroni               Show Bonferroni test line in p-value graphics.
    --fileformat=<type>        Formats available for the output: pdf, svg, jpg [default: pdf].
    --captions                 Generate figure captions.
    --outdir,-o=<dir>          Create a new directory if it doesn't exist. If not provided uses the input directory [default: graphics].

  """
    arg = docopt(mbsplot_doc, argv=None, help=True,
                 version='Maptools mbsplot version=0.1')
    #print(arg)
    df, arg = test_plot(arg)
    print(arg)
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

    if arg['--allele-freq-2'] == True:
        print('Allele Frequency 2 single chromosome')
        AF_mono_graph(df, arg, 'SNPidx2')

    if arg['--max-allele-freq2'] == True:
        AF_mono_graph(df, arg, 'MAX_SNPidx2')
        if arg['--multi-chrom'] == True:
           AF_multi_Vertical_graph(df, arg, 'MAX_SNPidx2')


def qtl_plot(argv):
    qtlplot_doc = """QTL-Plot
  Usage:
     maptools.py qtlplot <input_file> [options]
     maptools.py qtlplot -v
     maptools.py qtlplot -h

  Options:
    -h --help                  Show this screen.
    -v --version               Show the version.

  Output options:
    --chromosomes,-c=<opt>     List of chromosomes separeted by comma, from 1 to nºchr [default: all].
    --multi-chrom, -m          Multigraph representation.
    --pvalue, -p               Generates p-value graphics.
    --delta, -d                Generates the SNPidx for high and low pools and Delta graphic.
    --allele-freq-H, -H        Generates the alelle frequency graphics for the pool of high phenotype.
    --allele-freq-L, -L        Generates the allele frecuency graphics for the pool of low phenotype.
    --combine, -C              Combined multi-graphics with high and low lines. Use it together with --moving-avg.
    --euclidean-distance, -E   Generates graphics with Euclidean distance between individual SNPs and in groups of 100 adjacent markers.
    --g-statistic, -G          Generates graphics with G-statistic for individual SNPs.
    --qtl-seq, -Q              Generate muti-graphic with ED, G, DELTA and p-value graphics for each chromosome.
    --all, -a                  Generates all possible graphics.
    --moving-avg, -A=<int>     Selects the number of adjacent markers for the construction of moving average curves for allele frequncy and p-value.
    --palette, -P=<int>        Select the colour palette for your graphics (1 for colour blidness) [default: 1].
    --alpha, -t=<float>        Select the dots transparency on graphics (0.0 to 1.0) [default: 0.4].
    --bonferroni               Show Bonferroni test line in p-value graphics.
    --ci95                     Show the  95% cofindence interval in Delta graphics. 
    --fileformat=<type>        Formats available for the output: pdf, svg, jpg [default: pdf].
    --captions, -C             Generate figure captions.
    --outdir,-o=<dir>          Create a new directory if it doesn't exist. If not provided uses the input directory [default: graphics].

  """
    arg = docopt(qtlplot_doc, argv=None, help=True,version='qtl_plot qtl version=0.1')
    #print(arg)
    df, arg = test_plot(arg)
    print(arg)
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
   --mutagen, -M=<opt>           Type of mutagen used [defult: EMS].
   --region, -R=<region>         Region of the genome to explore (... -R chrName:Spos-Fpos)
Output Options:
   --output, -o=<file>           Output file.
   --outdir, -O=<dir>            Output directory [default: results].           
  """
  arg = docopt(annotate_doc, argv=None, help=True, version='Annotate variants version: 0.1')
  arg['pipe'] = sys.stdin.isatty()
  print(arg)
  arg = test_arg_ann(annotate_doc, arg)
  print(arg)
  df = create_df(arg)
  output = arg['--output']
  fsal = False
  arg['fsal'] = False
  if output != None:
    fsal = open(arg['--outdir']+arg['--output'], 'w')
    arg['fsal'] = fsal
  
  if arg['--input']:
    inp = open(arg['--input'], 'r')
  else:
    inp = sys.stdin
  for line in inp:
    if not line.startswith('#'):
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
  write_annotate_header(arg, fsal)
  for idx, row in df.iterrows():
    check_mutation2(row,arg)
    print('------------')
  finish = time.perf_counter()
  print(f'Variant annotation finished in {round(finish-start, 2)} second(s)')
  
  ##
