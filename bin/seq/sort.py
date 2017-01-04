#!/bin/env python
# Sort a SAM/BED/GTF/GFF file. Sorting is in-memory, so it is fast but inefficient.

import argparse
from collections import defaultdict, OrderedDict
import csv
from functools import cmp_to_key
import os
import sys
import re
from tqdm import tqdm
from xphyle import open_
from xphyle.utils import read_delimited, fileinput
from xphyle.paths import split_path

chr_map = dict(X=23, Y=24, M=25, Un=26)

def map_chrm(row, idx, chrm_map=None):
    if chrm_map is None:
        return row
    if row[idx] in chrm_map:
        row[idx] = chrm_map[row[idx]]
    else:
        return None

# Natural sorting of chromosomes, contigs, and other sequences
# that may be present in a GTF file
def compare_seq_names(a, b):
    def cmp(x, y):
        if x < y:
            return -1
        elif x == y:
            return 0
        else:
            return 1
    
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

class Seqname(object):
    def __init__(self, seqname):
        self.seqname = seqname
    
    def __lt__(self, other):
        return compare_seq_names(self.seqname, other.seqname) < 0
    
    def __repr__(self):
        return self.seqname

seqname_cache = {}
def seqname(name):
    if name not in seqname_cache:
        seqname_cache[name] = Seqname(name)
    return seqname_cache[name]

def lt(a, b, fields):
    for f in fields:
        a_val = getattr(a, f)
        b_val = getattr(b, f)
        if a_val < b_val:
            return True
        elif a_val > b_val:
            return False
    else:
        return False

class Record(object):
    def __init__(self, row, fields, cmp_fields):
        self.row = row
        self.fields = fields
        self.cmp_fields = cmp_fields
    
    def __lt__(self, other):
        return lt(self, other, self.cmp_fields)

    def __getattr__(self, name):
        field = self.fields[name]
        value = self.row[field[0]]
        type = field[1]
        if value == '.' and type != str:
            return None
        return type(value)

GFF_FIELDS = dict(
    seqname=(0, seqname),
    source=(1, str),
    feature=(2, str),
    start=(3, int),
    end=(4, int),
    score=(5, float),
    strand=(6, str),
    frame=(7, int),
    attribute=(8, str)
)
GFF_FEATURES = dict(
    gene=0,
    transcript=1,
    exon=2
)

class Gff(Record):
    def __init__(self, row):
        super().__init__(
            row, GFF_FIELDS, ('start', 'end', 'strand', 'feature_value'))
    
    def __getattr__(self, name):
        if name == 'feature_value':
            return GFF_FEATURES.get(self.feature, sys.maxsize)
        else:
            return super().__getattr__(name)

sequence_region_re = re.compile('##sequence-region ([^ ]+)')

def sort_gff(reader, output, chrm_map):
    seqs = defaultdict(lambda: [])
    file_headers = OrderedDict()
    chrm_headers = {}
    for line in reader:
        if line[0].startswith('#'):
            match = sequence_region_re.match(line[0])
            if match:
                chrm = match.group(1)
                if chrm_map and chrm in chrm_map:
                    chrm = chrm_map[chrm]
                chrm_headers[chrm] = line[0]
            else:
                file_headers[line[0]] = None
        else:
            line = map_chrm(line, 0, chrm_map)
            if not line:
                continue
            obj = Gff(line)
            seqs[obj.seqname].append(obj)
    
    # write the contents back to a file
    with open_(output, 'w') as o:
        w = csv.writer(o, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE,
                       doublequote=False)
        for line in file_headers.keys():
            w.writerow([line])
        for s in tqdm(sorted(seqs.keys()), desc="Sorting and writing rows...",
                      total=len(seqs)):
            if s.seqname in chrm_headers:
                w.writerow([chrm_headers[s.seqname]])
            rows = seqs[s]
            rows.sort()
            w.writerows(r.row for r in rows)

BED_FIELDS = dict(
    seqname=(0, seqname),
    start=(1, int),
    end=(2, int),
    name=(3, str),
    score=(4, float),
    strand=(5, str)
)

def sort_bed(reader, output, chrm_map):
    seqs = defaultdict(lambda: [])
    for line in reader:
        if line[0].startswith("#"):
            continue
        assert len(line) >= 3, "BED file must have at least 3 columns"
        line = map_chrm(line, 0, chrm_map)
        if not line:
            continue
        if len(line) < 6:
            line = line[0:3] + ['.','1000','.']
        obj = Record(line, BED_FIELDS, ('start', 'end', 'strand', 'name'))
        seqs[obj.seqname].append(obj)
    
    # write the contents back to a file
    with open_(output, 'w') as o:
        w = csv.writer(o, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE,
                       doublequote=False)
        for s in tqdm(sorted(seqs.keys()), desc="Sorting and writing rows...",
                      total=len(seqs)):
            rows = seqs[s]
            rows.sort()
            w.writerows(r.row for r in rows)

class Attr(object):
    def __init__(self, a):
        a = a.split(":")
        self.key = a[0]
        self.type = a[1]
        self.value = a[2]

SAM_FIELDS = dict(
    name=(0, str),
    flags=(1, int),
    seqname=(2, seqname),
    pos=(3, int),
    mapq=(4, int),
    cigar=(5, str),
    rnext=(6, str),
    pnext=(7, int),
    tlen=(8, int),
    seq=(9, str),
    qual=(10, str)
)

class Sam(Record):
    def __init__(self, row):
        super(Sam, self).__init__(row, SAM_FIELDS, ('seqname', 'pos', 'strand'))
        self._attributes = None
    
    @property
    def attributes(self):
        num_fields = len(self.row)
        if self._attributes is None and num_fields > 11:
            def parse_attr(a):
                attr = Attr(a)
                return (attr.key, attr)
            self._attributes = dict(
                parse_attr(self.row[i])
                for i in xrange(11, num_fields))
        return self._attributes
    
    @property
    def strand(self):
        a = self.attributes
        if a is not None and "XS" in a:
            return a["XS"].value
        return None

def sort_sam(reader, outfile, header=None, chrm_map=None):
    header_lines = list(read_delimited(header)) if header else None
    
    with open_(outfile, "w") as o:
        w = csv.writer(o, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE,
                       doublequote=False)
        if header_lines is not None:
            w.writerows(sort_sam_header(header_lines, chrm_map))
        
        seqs = defaultdict(lambda: [])
        for sam_line in reader:
            sam_line = map_chrm(sam_line, 2, chrm_map)
            if not sam_line:
                continue
            sam = Sam(row)
            seqs[sam.seqname].append(sam)
        
        for s in tqdm(sorted(seqs.keys()), desc="Sorting and writing rows...",
                      total=len(seqs)):
            rows = seqs[s]
            rows.sort(compare_sam_rows)
            w.writerows(r.row for r in rows)

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
        retval += sorted(seqname(x[1][3:]) for x in seq_headers)
    if "@RG" in h:
        retval += h["@RG"]
    if "@PG" in h:
        retval += h["@PG"]
    if "@CO" in h:
        retval += h["@CO"]
    return retval

def filter_sam(reader, outfile, chrm_list):
    """Only keep reads on chromosomes in chrm_list"""
    with open_(outfile, "w") as o:
        for sam_line in reader:
            w = csv.writer(o, delimiter="\t", quotechar="",
                           quoting=csv.QUOTE_NONE, doublequote=False)
            if sam_line[0].startswith("@"):
                if sam_line[0] == "@SQ" and sam_line[1][3:] not in chrm_list:
                    continue
                w.writerow(sam_line)
            else:
                sam = Sam(sam_line)
                if sam.seqname not in chrm_list:
                    continue
                w.writerow(sam.row)

FILE_TYPES = ("gff", "bed", "sam")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Performs proper sorting of BED, GFF, and SAM files")
    parser.add_argument("-T", "--type", choices=FILE_TYPES, default=None)
    parser.add_argument("-o", "--output", default="-",
        help='path name of properly sorted gtf')
    parser.add_argument("-H", "--sam_header", default=None,
        help='Header to prepend when sorting SAM files')
    parser.add_argument("-I", "--in_place", action="store_true", default=False,
        help="When a header file is specified, the output will be stored to "
             "that file")
    parser.add_argument("-m", "--chrm_map", default=None,
        help="Tab-delimited file with mapping of chromosome names in GTF/SAM "
             "to destination names")
    parser.add_argument("--filter_only", action="store_true", default=False,
        help="Only filter, don't sort.")
    parser.add_argument("infiles", nargs="*", default=['-'],
        help='paths of the files to sort')
    
    args = parser.parse_args()
    
    file_type = args.type
    if file_type is None:
        file_type = split_path(args.infiles[0], False)[-1]
    
    if file_type not in FILE_TYPES:
        raise Exception("Invalid file type: {0}".format(file_type))
    
    chrm_map = None
    if args.chrm_map is not None:
        l = list(read_delimited(args.chrm_map))
        if len(l[0]) == 1:
            chrm_map = dict((c[0], c[0]) for c in l)
        else:
            chrm_map = dict(l)
    
    reader = tqdm(read_delimited(fileinput(args.infiles)), desc="Reading rows...")
    
    if args.filter_only:
        filter_sam(reader, args.output, chrm_map)
    elif file_type == "gff":
        sort_gff(reader, args.output, chrm_map)
    elif file_type == "bed":
        sort_bed(reader, args.output, chrm_map)
    elif file_type == "sam":
        output = args.sam_header if args.sam_header and args.in_place else args.output
        sort_sam(reader, output, args.sam_header, chrm_map)
    else:
        raise Exception("Invalid file type {0}".format(args.type))
