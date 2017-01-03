#!/usr/bin/env python
# Subsample pairs of reads from a fastq file.
from argparse import ArgumentParser
import os
import random
from toolshed import nopen

def subsample(infiles, outfiles, prob, seed=None):
    prob = 1-prob

    if seed:
        random.seed(seed)

    def open_fq(f):
        fh = nopen(f, 'rb')
        return zip(*[fh] * 4)

    in_fh = [open_fq(i) for i in infiles]
    out_fh = [nopen(o, 'wb') for o in outfiles]

    try:
        written = 0
        for total, reads in enumerate(zip(*in_fh), 1):
            if random.random() >= prob:
                written += 1
                for read, fh in zip(reads, out_fh):
                    fh.writelines(read)

        print("wrote {} of {} reads".format(written, total))
    finally:
        for i in in_fh: i.close()
        for o in out_fh: o.close()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", default="{name}_{mate}.fastq.gz")
    parser.add_argument("-I", "--input-dir", default="")
    parser.add_argument("-o", "--output", default="{name}_{mate}.subsample.fastq.gz")
    parser.add_argument("-O", "--output-dir", default="")
    parser.add_argument("-m", "--mates", type=int, nargs=2, default=(1,2))
    parser.add_argument("-p", "--prob", default=0.1, type=float, help="Probability of retaining a read")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("name")
    args = parser.parse_args()

    infiles = [os.path.join(args.input_dir, args.input.format(name=args.name, mate=mate)) for mate in args.mates]
    outfiles = [os.path.join(args.output_dir, args.output.format(name=args.name, mate=mate)) for mate in args.mates]
    subsample(infiles, outfiles, args.prob, args.seed)
