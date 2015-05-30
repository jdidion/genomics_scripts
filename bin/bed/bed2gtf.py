#!/usr/bin/env python
# Converts a bed file to a valid gtf

import fileinput
import sys
from csv import reader, writer, QUOTE_NONE

if __name__ == "__main__":
    with open(sys.argv[2], 'w') as o:
        w = writer(o, delimiter="\t", quotechar="", quoting=QUOTE_NONE, doublequote=False)
        for i, row in enumerate(reader(fileinput.input(sys.argv[1]), delimiter="\t")):
            row = [row[0]] + ["bed2gtf", "exon"] + [int(row[1])+1] + [row[2]] + row[4:6] + [".", 'gene_id "{0}"'.format(row[3])]
            w.writerow(row)
