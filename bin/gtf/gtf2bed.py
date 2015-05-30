#!/usr/bin/env python
# Converts a gtf file to a valid bed

import fileinput
import sys
from csv import reader, writer, QUOTE_NONE

if __name__ == "__main__":
    with open(sys.argv[2], 'w') as o:
        w = writer(o, delimiter="\t", quotechar="", quoting=QUOTE_NONE, doublequote=False)
        for i, row in enumerate(reader(fileinput.input(sys.argv[1]), delimiter="\t")):
            row = [row[0], int(row[3]) - 1, row[4], "Locus{0}".format(i), 1000, row[6]]
            w.writerow(row)
