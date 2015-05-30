#!/usr/bin/env python
# For each sequence in a FASTA file, grep a FASTQ file (or other file with a set of sequences)
# for occurrences.

from Bio import SeqIO
import os
import sys

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
from util.cl import parse
from util.misc import bash

def main():
    def add_opts(parser):
        parser.add_argument("fasta", type="readable_file")
        parser.add_argument("seqs", type="readable_file")
        parser.add_argument("outdir", type="writeable_dir")
    
    args = parse(add_opts)
    
    with open(args.fasta, "rU") as f:
        for record in SeqIO.parse(f, "fasta"):
            outfile = os.path.join(args.outdir, "{0}.txt".format(record.id))
            cmd = 'grep "{0}" "{1}" > "{2}"'.format(record.seq, args.seqs, outfile)
            bash(cmd)

if __name__ == "__main__":
    main()
