#!/usr/bin/env python
# Read a FASTA file and output the C/G percentage for each entry

import os
import sys
from csv import writer

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
from util.cl import parse

def cg_pct(seq, outfile, window_size=100, decimal=False):
    with open(outfile, "w") as out:
        w = writer(out, delimiter="\t")
        seqlen = len(seq)
        for i in xrange(0, seqlen, window_size):
            end = min(i + window_size, seqlen)
            if i == end:
                break
            win = seq[i:end].upper()
            gc = float(sum(1 for char in win if char in ("C","G"))) / (end - i + 1)
            if not decimal:
                gc = int(round(gc * 100))
            w.writerow((i,end,gc))

def main(argv=None):
    def add_opts(parser):
        parser.add_argument("-d", "--decimal", action="store_true", default=False)
        parser.add_argument("-w", "--window_size", default=100, type=int)
        parser.add_argument("infile", type='readable_file')
        parser.add_argument("outfile", type='writeable_file')
    
    ns = parse(add_opts, args=argv)
    
    with open(ns.infile, 'rU') as inp:
        seq = inp.read().replace("\n","")

    cg_pct(seq, ns.outfile, ns.window_size, ns.decimal)

if __name__ == '__main__':
    main()