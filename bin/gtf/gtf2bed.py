#!/usr/bin/env python
# Converts a gtf file to a valid bed

from argparse import ArgumentParser
import fileinput
import re
import sys
from csv import reader, writer, QUOTE_NONE

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-n", "--name_column", type=int, default=None)
    parser.add_argument("-N", "--name_pattern", default=None)
    parser.add_argument("gtf")
    parser.add_argument("bed")
    args = parser.parse_args()
    
    name_column = args.name_column
    name_pattern = None if args.name_pattern is None else re.compile(args.name_pattern)
    
    with open(args.bed, 'w') as o:
        w = writer(o, delimiter="\t", quotechar="", quoting=QUOTE_NONE, doublequote=False)
        for i, row in enumerate(reader(fileinput.input(args.gtf), delimiter="\t")):
            
            name = "Locus{0}".format(i)
            if name_column is not None:
                name = row[args.name_column-1]
                if name_pattern is not None:
                    match = name_pattern.search(name)
                    if match is not None:
                        name = match.group(1)
                
            row = [row[0], int(row[3]) - 1, row[4], name, 1000, row[6]]
            w.writerow(row)
