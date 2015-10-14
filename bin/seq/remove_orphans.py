#!/usr/bin/env python

from argparse import ArgumentParser
from csv import reader, writer
import sys

def main():
    w = writer(sys.stdout, delimiter="\t")
    m1 = m2 = None

    for row in reader(sys.stdin, delimiter="\t"):
        if row[0][0] == "@":
            w.writerow(row)
        elif m1 is None:
            m1 = row
        elif m2 is None:
            if m1[0] == row[0]:
                m2 = row
            else:
                print "Removing orphaned read {0}".format(m1[0])
                m1 = row
        else:
            if m1[0] == row[0]:
                raise Exception("More than two reads with same ID")
            else:
                w.writerow(m1)
                w.writerow(m2)
                m1 = row
                m2 = None

    if m2 is not None:
        w.writerow(m1)
        w.writerow(m2)
    
if __name__ == "__main__":
    main()