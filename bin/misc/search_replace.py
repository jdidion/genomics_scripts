#!/usr/bin/env python
# Takes two files - one to have search and replace performed, the other with
# one line for each search/replace with the search term and replace term
# seaparated by a delimiter. An optional output file can be specified, otherwise
# the result will be written to the original file.

from csv import reader
import re
import os
import sys

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
from util.cl import parse
from util.io import csv_to_table, write_file

def main():
    def add_arguments(parser):
        parser.add_argument('-d', '--delim', default="\t")
        parser.add_argument('infile', type='readable_file')
        parser.add_argument('termfile', type='readable_file')
        parser.add_argument('outfile', type='writeable_file', nargs='?', default=None)
    ns = parse(add_arguments)
    
    with open(ns.infile, 'rU') as f: 
        s = f.read()
    
    for find, repl in csv_to_table(ns.termfile, ns.delim):
        s = re.sub(find, repl, s)
    
    outfile = ns.outfile or ns.infile
    write_file(outfile, s)
    
if __name__ == '__main__':
    main()
