#!/usr/bin/env python
# Identify polymorphic sites from an MSA in phylip format
import os
import sys

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
from util.io import safe_read_file_array

lines = safe_read_file_array(sys.argv[1])
names = [l[0:10] for l in lines[1:]]
seqs = [l[10:] for l in lines[1:]]
poss = []

for pos in xrange(0, len(seqs[0])):
    bases = [s[pos] for s in seqs]
    s = set(bases)
    if len(s) > 1 and 'N' not in s:# (len(s) > 2 or 'N' not in s):
        poss.append(pos)

print "%i %i" % (len(seqs), len(poss))
for l in xrange(0, len(seqs)):
    print names[l] + ''.join(seqs[l][p] for p in poss)
    