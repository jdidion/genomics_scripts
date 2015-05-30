#!/bin/env python
# Sort a SAM or GTF file. Sorting is in-memory, so it is fast but inefficient.

import os
import sys
import csv
import re
import argparse
import fileinput
from itertools import tee, imap
from collections import defaultdict, namedtuple

from jpd.util.io import open_mgr

chr_map = dict(X=23, Y=24, M=25, Un=26)

# Natural sorting of chromosomes, contigs, and other sequences
# that may be present in a GTF file
def compare_seq_names(a, b):
    if a == b:
        return 0
    
    # strip off chr prefix
    if a.startswith("chr"):
        a = a[3:]
    if b.startswith("chr"):
        b = b[3:]
    
    # map sex chromosome names to numbers (human genome specific!)
    a = chr_map.get(a, a)
    b = chr_map.get(b, b)
    a_int = b_int = False
    
    # try to convert to integers
    try:
        a = int(a)
        a_int = True
    except:
        pass
    
    try:
        b = int(b)
        b_int = True
    except:
        pass
    
    # natural sorting of numbers, with numbered chromosomes
    # always coming before contigs, etc.
    if a_int and b_int:
        return cmp(a, b)
    elif a_int:
        return -1
    elif b_int:
        return 1
    else:
        # if these are contigs that have been localized to a chromosome,
        # sort by the chromosome
        ai = a.find("_")
        bi = b.find("_")
        if ai >= 0 and bi >= 0:
            a_chr = a[0:ai]
            b_chr = b[0:bi]
            c = compare_seq_names(a_chr, b_chr)
            if c != 0:
                return c
        elif ai >= 0:
            return -1
        elif bi >= 0:
            return 1
    
    # default to string sorting
    return cmp(a, b)

# compare two rows on the same chromosome
def compare_rows(a, b):
    c = cmp(a.start, b.start)
    if c == 0:
        c = cmp(a.end, b.end)
    if c == 0:
        c = cmp(a.strand, b.strand)
    return c

def sort(infiles, output, chrm_map, make_fn, cmp_fn):
    seqs = defaultdict(lambda: [])
    for line in csv.reader(fileinput.input(infiles, mode="rU"), delimiter='\t'):
        if line[0].startswith("#"):
            continue
        obj = make_fn(line, chrm_map)
        if obj is not None:
            seqs[obj.seqname].append(obj)
    
    # write the contents back to a file
    with open_mgr(output, 'wb') as o:
        w = csv.writer(o, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE, doublequote=False)
        for s in sorted(seqs.keys(), compare_seq_names): 
            rows = seqs[s]
            rows.sort(cmp_fn)
            w.writerows(rows)

Gtf = namedtuple("Gtf", 
    ("seqname","source","feature","start","end","score","strand","frame","attribute"))
def make_gtf(row, chrm_map=None):
    if chrm_map is not None:
        if row[0] in chrm_map:
            row[0] = chrm_map[row[0]]
        else:
            return None
    row[3] = int(row[3])
    row[4] = int(row[4])
    return(Gtf(*row))

GTF_FEATURES = dict(gene=0, transcript=1, exon=2)
def compare_gtf_rows(a, b):
    c = compare_rows(a, b)
    if c == 0:
        af = GTF_FEATURES.get(a.feature, 3)
        bf = GTF_FEATURES.get(b.feature, 3)
        c = cmp(af, bf)
    return c

def sort_gtf(infiles, output, chrm_map=None):
    sort(infiles, output, chrm_map, make_fn=make_gtf, cmp_fn=compare_gtf_rows)

Bed = namedtuple("Bed", 
    ("seqname","start","end","name","score","strand"))
def make_bed(row, chrm_map=None):
    assert len(row) >= 3, "BED file must have at least 3 columns"
    if chrm_map is not None:
        if row[0] in chrm_map:
            row[0] = chrm_map[row[0]]
        else:
            return None
    row[1] = int(row[1])
    row[2] = int(row[2])
    if len(row) < 6:
        row = row[0:3] + ['.','1000','.']
    elif len(row) > 6:
        row = row[0:6]
    return(Bed(*row))

def compare_bed_rows(a, b):
    c = compare_rows(a, b)
    if c == 0:
        c = cmp(a.name, b.name)
    return c

def sort_bed(infiles, output, chrm_map=None):
    sort(infiles, output, chrm_map, make_fn=make_bed, cmp_fn=compare_bed_rows)

class Attr(object):
    def __init__(self, a):
        a = a.split(":")
        self.key = a[0]
        self.type = a[1]
        self.value = a[2]

def make_sam(row, chrm_map=None):
    if chrm_map is not None:
        if row[2] in chrm_map:
            row[2] = chrm_map[row[2]]
        else:
            return None
    return Sam(row)

class Sam(object):
    def __init__(self, row):
        self.row = row
        self._attributes = None
        
    @property
    def seqname(self): return self.row[2]
        
    @property
    def pos(self): return int(self.row[3])
    
    @property
    def attributes(self):
        N = len(self.row)
        if self._attributes is None and N > 11:
            def parse_attr(a):
                attr = Attr(a)
                return (attr.key, attr)
            self._attributes = dict(parse_attr(self.row[i]) for i in xrange(11, N))
        return self._attributes
    
    @property
    def strand(self):
        a = self.attributes
        if a is not None and "XS" in a:
            return a["XS"].value
        return None

def sort_sam_header(rows, chrm_map=None):
    h = defaultdict(lambda: [])
    for r in rows:
        h[r[0]].append(r)
    retval = []
    if "@HD" in h:
        retval += h["@HD"]
    if "@SQ" in h:
        seq_headers = h["@SQ"]
        seq_names = list(x[1][3:] for x in seq_headers)
        if chrm_map is not None:
            new_seq_headers = []
            for i,sn in enumerate(seq_names):
                if sn in chrm_map:
                    sh = seq_headers[i]
                    sh[1] = "SN:{0}".format(chrm_map[sn])
                    new_seq_headers.append(sh)
            seq_headers = new_seq_headers
        retval += sorted(seq_headers, compare_seq_names, lambda x: x[1][3:])
    if "@RG" in h:
        retval += h["@RG"]
    if "@PG" in h:
        retval += h["@PG"]
    if "@CO" in h:
        retval += h["@CO"]
    return retval

def compare_sam_rows(a, b):
    c = cmp(a.pos, b.pos)
    if c == 0:
        c = cmp(a.strand, b.strand)
    return c

def sort_sam(infiles, outfile, header=None, chrm_map=None):
    header_lines = None
    if header is not None:
        with open(header, "rU") as i:
            header_lines = list(csv.reader(i, delimiter='\t'))
    
    with open_mgr(outfile, "w") as o:
        w = csv.writer(o, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE, doublequote=False)
        if header_lines is not None:
            w.writerows(sort_sam_header(header_lines, chrm_map))
        
        seqs = defaultdict(lambda: [])
        for sam_line in csv.reader(fileinput.input(infiles), delimiter='\t'):
            sam = make_sam(sam_line, chrm_map)
            if sam is not None:
                seqs[sam.seqname].append(sam)
        
        for s in sorted(seqs.keys(), compare_seq_names): 
            rows = seqs[s]
            rows.sort(compare_sam_rows)
            for r in rows:
                w.writerow(r.row)

def filter_sam(infiles, outfile, chrm_list):
    """Only keep reads on chromosomes in chrm_list"""
    with open_mgr(outfile, "w") as o:
        for sam_line in csv.reader(fileinput.input(infiles), delimiter='\t'):
            w = csv.writer(o, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE, doublequote=False)
            if sam_line[0].startswith("@"):
                if sam_line[0] == "@SQ" and sam_line[1][3:] not in chrm_list:
                    continue
                w.writerow(sam_line)
            else:
                sam = Sam(sam_line)
                if sam.seqname not in chrm_list:
                    continue
                w.writerow(sam.row)

FILE_TYPES = ("gtf","bed","sam")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Performs proper sorting of GTF""")
    parser.add_argument("-T", "--type", choices=FILE_TYPES, default=None)
    parser.add_argument("-i", "--input", default="-", help='path to gtf file to sort')
    parser.add_argument("-o", "--output", default="-", help='path name of properly sorted gtf')
    parser.add_argument("-H", "--sam_header", default=None, help='Header to prepend when sorting SAM files')
    parser.add_argument("-I", "--in_place", action="store_true", default=False, 
        help="When a header file is specified, the output will be stored to that file")
    parser.add_argument("-m", "--chrm_map", default=None,
        help="Tab-delimited file with mapping of chromosome names in GTF/SAM to destination names")
    parser.add_argument("--filter_only", action="store_true", default=False,
        help="Only filter, don't sort.")
    parser.add_argument("infiles", nargs="*", help='path to gtf file to sort')
    
    args = parser.parse_args()
    
    file_type = args.type
    if file_type is None:
        file_type = os.path.splitext(args.infiles[0])[1][1:]
    
    if file_type not in FILE_TYPES:
        raise Exception("Invalid file type: {0}".format(file_type))    
    
    chrm_map = None
    if args.chrm_map is not None:
        with open(args.chrm_map, "rU") as m:
            l = list(csv.reader(m, delimiter="\t"))
            if len(l[0]) == 1:
                chrm_map = dict((c[0],c[0]) for c in l)
            else:
                chrm_map = dict(l)
    
    infiles = args.infiles if len(args.infiles) > 0 else [args.input]
    
    if args.filter_only:
        filter_sam(infiles, args.output, chrm_map)
    elif file_type == "gtf":
        sort_gtf(infiles, args.output, chrm_map)
    elif file_type == "bed":
        sort_bed(infiles, args.output, chrm_map)
    elif file_type == "sam":
        output = args.sam_header if args.sam_header and args.in_place else args.output
        sort_sam(infiles, output, args.sam_header, chrm_map)
    else:
        raise Exception("Invalid file type {0}".format(args.type))