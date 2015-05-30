#!/usr/bin/env python
# Convert a FASTA file into CSV format (header,sequence).

from Bio import SeqIO
import os
import re
import sys

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
from util.cl import parse

def convert(fasta, out, delim, replace):
    for entry in SeqIO.parse(fasta, "fasta"):
        header = re.sub("[%s]" % replace, delim, entry.name)
        seq = entry.seq.upper()
        out.write(delim.join((header, str(seq))))
        out.write("\n")

def main(argv=None):
    def add_opts(parser):
        parser.add_argument("-d", "--delim", default=",")
        parser.add_argument("-r", "--replace", default=":-")
        parser.add_argument("infile", type='readable_file')
        parser.add_argument("outfile", type='writeable_file')
    
    ns = parse(add_opts, args=argv)
    
    with open(ns.infile, 'rU') as i, open(ns.outfile, 'w') as o:
        convert(i, o, ns.delim, ns.replace)

if __name__ == '__main__':
    main()