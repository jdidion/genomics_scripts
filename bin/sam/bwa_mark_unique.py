#!/bin/bash
# Set a flag on every read in a bam file.

from argparse import ArgumentParser
import pysam

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-t", "--threshold", type=int, default=1,
        help="Minimum difference between AS and XS for the read to be marked unique.")
    parser.add_argument("infile")
    parser.add_argument("outfile")
    args = parser.parse_args()
    
    infile = pysam.AlignmentFile(args.infile, "rb")
    outfile = pysam.AlignmentFile(args.outfile, "wb", template=infile)
    
    for read in infile:
        xt = 'N'
        if read.has_tag('AS'):
            if read.has_tag('XS'):
                if read.get_tag('AS') - read.get_tag('XS') >= args.threshold:
                    xt = 'U'
                else:
                    xt = 'R'
            else:
                xt = 'U'
        read.set_tag('XT', xt, 'Z')
        outfile.write(read)
    
    infile.close()
    outfile.close()
