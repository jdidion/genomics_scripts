#!/usr/bin/env python
# Extracts IDs from cuffnorm summary count/fpkm tables based on minimum criteria.

from csv import reader
import sys

# TODO: add filtering options; right now it includes any row with at least one column with FPKM >= 1

def main():
    with open(sys.argv[1], 'rU') as i, open(sys.argv[2], 'w') as o:
        r = reader(i, delimiter='\t')
        r.next() # skip header
        for row in r:
            for col in row[1:]:
                if float(col) >= 1:
                    o.write(row[0])
                    o.write('\n')
        

if __name__ == "__main__":
    main()