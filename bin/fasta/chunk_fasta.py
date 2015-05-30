#!/usr/bin/env python
# Split records in a FASTA file into smaller chunks of a fixed size.

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys

def chunk_iter(i, size):
    for record in SeqIO.parse(handle, "fasta"):
        rec_id = record.id
        seq_len = len(record.seq)
        seq_type = record.seq.alphabet
        for idx, start in enumerate(xrange(0, len(record.seq), size)):
            end = min(start+size, seq_len)
            yield SeqRecord(Seq(record.seq[start:end], seq_type), 
                id="{0}_{1}".format(rec_id, idx))

def chunk(i, o, size=100):
    SeqIO.write(chunk_iter(i, size), o, "fasta")
        
def main():
    parser = ArgumentParser()
    parser.add_argument("-s", "--size", type=int, default=100)
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()

    with open(args.input, "rU") as i, open(args.output, "w") as o:
        chunk(i, o, args.size)

if __name__ == "__main__":
    main()