#!/usr/bin/env python
# Extracts rows from a GTF that match a list of IDs.

from csv import reader
import sys
import re

# TODO: add filtering options

def main():
    with open(sys.argv[1], 'rU') as i:
        ids = set(i)
    
    pat = re.compile('(TCONS_\d+)')
    with open(sys.argv[2], 'rU') as i, open(sys.argv[3], 'w') as o:
        for line in i:
            match = pat.search(line)
            if match is not None and match.group(1) in ids:
                o.write(line)

if __name__ == "__main__":
    main()