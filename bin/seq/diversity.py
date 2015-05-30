#!/usr/bin/env python
# Calculate sequence diversity metrics for a set of sequences. Currently, this only handles a
# stream of sequences (one per line) of equal size. This can easily be extracted from a BAM
# file using samtools: samtools foo.bam | cut -f | diversity.py <options> -. For
# every kmer size, one summary file and one counts file is written (unless --summary_only is
# specified) to the output directory. If you will be processing tens-of-billions of reads or
# more, you should specify the --bigmem argument.

import fileinput
import os
import sys
import numpy as np

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
from util.cl import parse
from util.gnet import SeqDiv

def main():
    def add_args(parser):
        parser.add_argument("-k", "--kmer_sizes", type="int_list", action="extend_overwrite", default=(1,))
        parser.add_argument("-r", "--read_length", type=int, default=100)
        parser.add_argument("--bigint", action="store_true", default=False)
        parser.add_argument("--log", type="writeable_file", default=None)
        parser.add_argument("--log_interval", type=int, default=100000)
        parser.add_argument("--prefix", default=None)
        parser.add_argument("--summary_only", action="store_true", default=False)
        parser.add_argument("input_file", type="readable_file")
        parser.add_argument("output_dir", type="writeable_dir", nargs="?", default=".")
    
    args = parse(add_args)
    
    dtype = np.int64 if args.bigint else np.int32
    
    total_mem = sum(SeqDiv.estimate_memory(args.read_length, k, dtype) for k in args.kmer_sizes) / 1000000.0
    sys.stderr.write("This program will use up to {0} MiB of memory\n".format(total_mem))
    
    kmers = dict((k, SeqDiv(k, dtype)) for k in args.kmer_sizes)
    
    log = None
    start = None
    if args.log is not None:
        from datetime import datetime
        start = datetime.now()
        log = open(args.log, "w", 0)
        log.write("Starting at {0}\n".format(start))
        log.flush()
        
    for read_num, read in enumerate(fileinput.input(args.input_file, mode="rU"), 1):
        # TODO: this could be threaded if it's too slow
        read = read.strip()
        read_len = len(read)
        if read_len != args.read_length:
            sys.exit("Invalid read length at read {0}: expected {1}, actual {2}".format(
                read_num, args.read_length, read_len))
        for sd in kmers.values():
            sd.insert_nocheck(read, read_len)
        
        if log is not None and read_num % args.log_interval == 0:
            now = datetime.now()
            log.write("Processed {0} reads in {1} hours\n".format(read_num, 
                round((now - start).total_seconds() / 3600, 3)))
            log.flush()
    
    prefix = args.prefix
    if prefix is None:
        prefix = os.path.splitext(os.path.basename(args.input_file))[0]
    
    for k in kmers.keys():
        summary_file = os.path.join(args.output_dir, "{0}_summary_{1}mers.csv".format(prefix, k))
        with open(summary_file, "w") as o:
            kmers[k].write_summary(o)
        if not args.summary_only:
            count_file = os.path.join(args.output_dir, "{0}_counts_{1}mers.csv".format(prefix, k))
            with open(count_file, "w") as o:
                kmers[k].write_counts(o)
    
    if log is not None:
        log.close()

if __name__ == "__main__":
    main()