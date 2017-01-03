#!/usr/bin/env python
# Split fastq into N roughly equal-sized files
# Note: this requires py35 (otherwise substitute izip for zip)
from argparse import ArgumentParser
import gzip
import os
import random
from toolshed import nopen

def subsample(infile, outfiles):
    in_fh = nopen(infile, 'rb')
    out_fh = [gzip.open(o, 'wt') for o in outfiles]
    n = len(out_fh)

    try:
        for i, read in enumerate(zip(*[in_fh] * 4)):
            out_fh[i % n].writelines(read)
    finally:
        in_fh.close()
        for o in out_fh: o.close()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-O", "--output-dir", default=None)
    parser.add_argument("-p", "--parts", type=int, default=2)
    parser.add_argument("infile")
    args = parser.parse_args()
    
    inpath = os.path.abspath(args.infile)
    infile = os.path.basename(inpath)
    output_dir = args.output_dir or os.path.dirname(args.infile)
    outfiles = [
        os.path.join(output_dir, "{}.{}".format(infile, index))
        for index in range(args.parts)
    ]
    subsample(args.infile, outfiles)
