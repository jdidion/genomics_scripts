#!/usr/bin/env python

from csv import reader, writer
import sys

def main():
    w = writer(sys.stdout, delimiter="\t")
    for i,row in enumerate(reader(sys.stdin, delimiter="\t")):
        if row[0] == 'spliceName':
            continue
        if len(row) == 6:
            w.writerow((row[1], int(row[3]) - 1, row[4], row[0], row[5], row[2]))
        else:
            w.writerow((row[0], int(row[2]) - 1, row[3], 'Novel{0}'.format(i), row[4], row[1]))

if __name__ == "__main__":
    main()