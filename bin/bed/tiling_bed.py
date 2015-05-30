#!/usr/bin/env python
# Create a BED file that tiles the genome with windows of a specified size. The Y chromosome
# won't be tiled, and a separate window size can be specified for chrM.

import fileinput
import os
import sys
import numpy as np

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
from util.cl import parse
from util.io import csv_to_dict

def main():
    def add_args(parser):
        parser.add_argument("-c", "--chromosomes", action="extend_overwrite", type="str_list", 
            default=[str(c) for c in xrange(1,20)] + ["X","M"])
        parser.add_argument("-m", "--mt_window_size", type=int, default=100)
        parser.add_argument("-w", "--window_size", type=int, default=1000000)
        parser.add_argument("--start", type=int, default=3000000,
            help="Start of first window for autosomal and X chromosomes.")
        parser.add_argument("genome_file", type="readable_file")
        parser.add_argument("output_file", type="writeable_file")
    
    args = parse(add_args)
    
    genome = csv_to_dict(args.genome_file, delim="\t")
    
    with open(args.output_file, "w") as o:
        for chrm in args.chromosomes:
            chrm_name = "chr{0}".format(chrm)
            chrm_size = int(genome[chrm_name][1])
            window_size = args.mt_window_size if chrm == "M" else args.window_size
            
            for pos in xrange(args.start, chrm_size, window_size):
                if pos < chrm_size:
                    o.write("{0}\t{1}\t{2}\n".format(chrm_name, pos, min(pos + window_size, chrm_size)))
        
if __name__ == "__main__":
    main()