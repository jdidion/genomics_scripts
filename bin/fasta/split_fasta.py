#!/usr/bin/env python

import os
import sys
import fileinput

sys.path.append(os.path.join(os.environ['LAB_HOME'], 'lib/python'))
from util.cl import parse


def main():
    def add_args(parser):
        parser.add_argument('infile', type='readable_file', metavar='FILE')
        parser.add_argument('outdir', type='writeable_dir', default='.', nargs='?')
    
    ns = parse(add_args)
    
    outfile = None
    for line in fileinput.input(ns.infile):
        if line[0] == ">":
            if outfile is not None:
                outfile.close()
            fname = os.path.join(ns.outdir, "{0}.fasta".format(line[1:-1]))
            outfile = open(fname, 'w')
        outfile.write(line)

if __name__ == '__main__':
    main()