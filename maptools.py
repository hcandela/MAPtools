#!/usr/bin/env python3
"""
Usage:
   maptools.py <command> <args>...
   maptools.py <command>
   maptools.py

Available commands:
   mbs
   mbsplot
   qtl
   qtlplot
   merge
   annotate

Version:
   Maptools version: 1.00
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