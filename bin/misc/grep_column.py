#!/usr/bin/env python

from argparse import ArgumentParser
from csv import reader, writer

def main():
    parser = ArgumentParser()
    parser.add_argument('-c', '--column', type=int, default=1)
    parser.add_argument('-f', '--file', help="File with strings to match, one per line")
    parser.add_argument('-o', '--output', help="Output file (defaults to stdout)")
    parser.add_argument('input', help="File to search")
    args = parser.parse_args()
    
    strs = set(s.rstrip() for s in open(args.file, 'rU'))
    with open(args.input, 'rU') as i, open(args.output, 'w') as o:
        w = writer(o, delimiter="\t")
        for row in reader(i, delimiter="\t"):
            if row[args.column-1] in strs:
                w.writerow(row)

if __name__ == "__main__":
    main()