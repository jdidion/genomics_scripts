#!/usr/bin/env python

from collections import defaultdict
from csv import writer
import pysam
import sys

sam = pysam.Samfile(sys.argv[1], "rb")
hist = defaultdict(lambda:0)

pos = None
cnt = 0
for read in sam.fetch():
    if pos == read.pos:
        cnt += 1
    else:
        if pos is not None:
            hist[cnt] += 1
        pos = read.pos
        cnt = 1

with open(sys.argv[2], "w") as o:
    w = writer(o)
    for k,v in hist.iteritems():
        w.writerow((v,k))