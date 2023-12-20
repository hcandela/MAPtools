#!/usr/bin/env python3
"""
Program: MAPtools
Version: 1.00
Contact: hcandela@umh.es

Usage:
   maptools.py <command> <args>...
   maptools.py <command>
   maptools.py

Available commands:
   mbs
   qtl
   merge
   plot
   annotate
   citation
"""

from docopt import docopt
from pylibs.commands import *

if __name__ == '__main__':
    args = docopt(__doc__, version=v_maptools, options_first=True)
    if not args['<command>']:
        print(__doc__,end='', file=sys.stderr)
        sys.exit()
    if not args['<args>']:
        argv = [args['<command>']]
    else:
        argv = [args['<command>']] + args['<args>']
    if args['<command>'] == 'mbs':
        mbs(argv)
    elif args['<command>'] == 'qtl':
        qtl(argv)
    elif args['<command>'] == 'plot':
        plot(argv)
    elif args['<command>'] == 'merge':
        merge(argv)
    #elif args['<command>'] == 'annotatevcf':
    #    annotatevcf(argv)
    elif args['<command>'] == 'annotate':
        annotate(argv)
    elif args['<command>'] == 'citation':
        citation(argv)
    else:
        print('%r is not a maptools command. See \'maptools --help\'.'% args['<command>'], file=sys.stderr)
        sys.exit()