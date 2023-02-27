#!/usr/bin/env python3
""" MBS & QTL-seq Analysis
About:
   Variant analyzer for MBS & QTL-seq experiment from VCF files. To be used in conjunction with bcftools mpileup and
   bcftools call. Also it can be used giving the vcf input on stream with cat.
   First option:
   bcftools mpileup -f ref.fa --annotate FORMAT/AD b1.bam b2.bam ... | bcftools call (-m | -mv) -V indels | python3 maptools.py <command> <args>...
   Second option:
   cat mpileup_call_in.vcf | python3 maptools.py <command> <args>...
   Third option:
   python3 maptools.py <command> <mpileup_call_in.vcf> <args>...
Usage:
   maptools.py <command> <args>...
   maptools.py <command>
   maptools.py --version
   maptools.py -h
   maptools.py

Availaible commands:
   mbs
   mbsplot
   qtl
   qtlplot
   merge
   annotate
"""

from docopt import docopt
from pylibs.commands import *

if __name__ == '__main__':
    args = docopt(__doc__, version='maptools version=0.1', options_first=True)
    if not args['<command>']:
        print(__doc__,end='')
        sys.exit()
    if not args['<args>']:
        argv = [args['<command>']]
    else:
        argv = [args['<command>']] + args['<args>']
    if args['<command>'] == 'mbs':
        mbs(argv)
    elif args['<command>'] == 'qtl':
        qtl(argv)
    elif args['<command>'] == 'mbsplot':
        mbs_plot(argv)
    elif args['<command>'] == 'qtlplot':
        qtl_plot(argv)
    elif args['<command>'] == 'merge':
        merge(argv)
    elif args['<command>'] == 'annotate':
        annotate(argv)
    else:
        sys.exit('%r is not a maptools command. See \'maptools --help\'.'% args['<command>'])