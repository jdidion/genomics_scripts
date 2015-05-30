#!/usr/bin/env python
"""
This program takes a file containing fragments and another containing MSRE sites and
writes to a third file a line for each MSRE site found within each fragment. It is assumed
that both input files are sorted and that there is no possible overlap of cut sites with
recognition sequences.

The output file looks like::

    proximal_NspI,distal_NspI,MSRE_site
"""

import os
import sys
import csv

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
import util.cl

def int_file_generator(f):
    for line in f:
        yield int(line.strip())

def find_cuts_in_fragments(frag_file, msre_file, output_file):
    with open(frag_file) as frags, open(msre_file) as msre, open(output_file, 'w') as out:
        frag_reader = csv.reader(frags)
        msre_reader = int_file_generator(msre)
        
        cur_site = None
        for frag in frag_reader:
            start, end = map(int, frag)

            while cur_site < start:
                try:
                    cur_site = msre_reader.next()
                except:
                    return
                        
            if start < cur_site < end:
                out.write("%i,%i,%i\n" % (start,end,cur_site))

if __name__ == '__main__':
    def add_args(parser):
        parser.add_argument("frag_file", metavar="FILE", type="readable_file")
        parser.add_argument("msre_file", metavar="FILE", type="readable_file")
        parser.add_argument("output_file", metavar="FILE", type="writeable_file")
    
    ns = util.cl.parse(add_args)
    find_cuts_in_fragments(ns.frag_file, ns.msre_file, ns.output_file)
    