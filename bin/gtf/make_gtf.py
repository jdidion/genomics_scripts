#!/usr/bin/env python
# Converts a simple chrm,start,end,strand table to a valid bedfile.

import fileinput
import sys
from csv import reader, writer

if __name__ == "__main__":
    with open(sys.argv[2], 'w') as o:
        w = writer(o, delimiter="\t")
        for i, row in enumerate(reader(fileinput.input(sys.argv[1]), delimiter="\t")):
            row = [row[0], "make_gtf", "exon", int(row[1]) + 1, row[2], ".", row[3], ".", row[4]]
            w.writerow(row)
