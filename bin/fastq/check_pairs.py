#!/usr/bin/env python
# Given a pair of fastq files, ensure that each has the exact same set of reads
# in the same order. Optionally, estimate genome coverage from the total number
# of reads * read length / genome size.

from Bio import SeqIO
import argparse
import csv
import gzip
import os
import re
import time
progress = True
try:
    import progressbar
    import progressbar.widgets as widgets
except:
    progress = False

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-z", "--gzip", action="store_true", default=False,
        help="Files are gzipped, even if they don't have a .gz extension")
    parser.add_argument("-g", "--genome-size", default=3137161264) # from hg19
    parser.add_argument("-p", "--name-pattern", default=".* (.*) length=(\d+)")
    parser.add_argument("-o", "--output-file", default=None, help="Write summary to tsv file "\
        "(name, min_read_length, max_read_length, read_count, total_Gb, coverage_est)")
    parser.add_argument("-q", "--quiet", action="store_true", default=False,
        help="Don't print summary to stdout")
    parser.add_argument("-s", "--sample_name", default=None,
        help="Sample name; otherwise file name will be used")
    parser.add_argument("fastq1")
    parser.add_argument("fastq2")
    args = parser.parse_args()

sample = args.sample_name or re.match("(.*)_?1.f(ast)?q(.gz)?", args.fastq1).group(1)
name_pattern = re.compile(args.name_pattern)

def parse_desc(r):
    desc = name_pattern.match(r.description)
    if desc is None:
        raise Exception("Name pattern {} does not match read description {}".format(
            args.name_pattern, r.description))
    g = desc.groups()
    if len(g) == 1:
        name = g[0]
        length = len(r.seq)
    else:
        name, length = g
        length = int(length)
    return (name, length)

gz = args.gzip or args.fastq1.endswith(".gz")
if gz:
    fq1 = gzip.open(args.fastq1, "rt")
    fq2 = gzip.open(args.fastq2, "rt")
else:
    fq1 = open(args.fastq1, "rU")
    fq2 = open(args.fastq2, "rU")

read_length = None
min_read_length = None
max_read_length = None
total_bp = None
read_count = 0

print("Reading from paired-end files for sample {}...".format(sample))

start = time.process_time()

with fq1, fq2:
    itr = zip(SeqIO.parse(fq1, "fastq"), SeqIO.parse(fq2, "fastq"))
    if progress:
        pb = progressbar.ProgressBar(widgets = [
            widgets.AnimatedMarker(),
            ' ', widgets.Counter(),
            ' read pairs ; ', widgets.Timer(),
        ], poll_interval=1000)
        itr = pb(itr)
    for read1, read2 in itr:
        name1, length1 = parse_desc(read1)
        name2, length2 = parse_desc(read2)
        if name1 != name2:
            raise Exception("Unpaired reads: {0} <> {1}".format(name1, name2))
        cur_read_length = length1 + length2
        
        if total_bp is None and read_length is not None and read_length != cur_read_length:
            print("Read lengths are variable: {} != {}".format(read_length, cur_read_length))
            if read_length is None:
                total_bp = 0
            else:
                total_bp = read_length * read_count
                read_length = None
        
        read_count += 1
        
        if total_bp is not None:
            total_bp += cur_read_length
            min_read_length = cur_read_length if min_read_length is None else min(min_read_length, cur_read_length)
            max_read_length = cur_read_length if max_read_length is None else max(max_read_length, cur_read_length)
        
        elif read_length is None:
            min_read_length = max_read_length = read_length = cur_read_length

elapsed_time = time.process_time() - start

if total_bp is None:
    total_bp = read_length * read_count

total_gb = total_bp / 1000000000.0
coverage = total_gb / args.genome_size

if args.output_file is not None:
    with open(args.output_file, "w") as o:
        csv.writer(o).writerow((sample, min_read_length, max_read_length, read_count, total_gb, coverage))
    
if not args.quiet:
    print("Read length: {}".format("{}-{}".format(min_read_length, max_read_length) 
        if min_read_length != max_read_length else read_length))
    print("Total read count: {}".format(read_count))
    print("Total base count (in Gb): {}".format(round(total_gb, 2)))
    print("Approx genome coverage: {}X".format(round(coverage, 1)))
    print("Elapsed time: {}s".format(round(elapsed_time, 1)))
