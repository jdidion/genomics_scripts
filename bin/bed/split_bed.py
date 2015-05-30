#!/usr/bin/env python
# Splits a bed file with block definitions into one row per block

import argparse
import fileinput
import sys
from csv import reader, writer

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("outfile")
    args = parser.parse_args()
    
    with open(args.outfile, 'w') as o:
        w = writer(o, delimiter="\t")
        for i, row in enumerate(reader(fileinput.input(args.infile), delimiter="\t")):
            loc_start = int(row[1])
            sizes = map(int, row[10].split(','))
            starts = map(int, row[11].split(','))
            for i in xrange(len(sizes)):
                block = [row[0], loc_start + starts[i], loc_start + starts[i] + sizes[i]] + row[3:6]
                w.writerow(block)
