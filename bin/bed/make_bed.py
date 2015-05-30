#!/usr/bin/env python
# Converts a simple chrm,start,end,strand table to a valid bedfile.

import argparse
import fileinput
import sys
from csv import reader, writer

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--one_indexed", action="store_true", default=False)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("infile")
    parser.add_argument("outfile")
    args = parser.parse_args()
    
    with open(args.outfile, 'w') as o:
        w = writer(o, delimiter="\t")
        for i, row in enumerate(reader(fileinput.input(args.infile), delimiter="\t")):
            if i < args.skip:
                continue
            if args.one_indexed:
                row = [row[0], int(row[1])-1, row[2], "Locus{0}".format(i), 1000, row[3]]
            else:
                row = row[0:3] + ["Locus{0}".format(i), 1000] + [row[3]]
            w.writerow(row)
