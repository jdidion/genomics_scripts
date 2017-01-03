#!/bin/bash
# Set a flag on every read in a bam file.

from argparse import ArgumentParser
import pysam

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--flag", type=int, default=1)
    parser.add_argument("infile")
    parser.add_argument("outfile")
    args = parser.parse_args()
    
    infile = pysam.AlignmentFile(args.infile, "rb")
    outfile = pysam.AlignmentFile(args.outfile, "wb", template=infile)
    
    for read in infile:
        read.flag |= args.flag
        outfile.write(read)
    
    infile.close()
    outfile.close()
