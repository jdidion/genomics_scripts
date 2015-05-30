#!/usr/bin/env python
# Show the genotypes for a given variant

import vcf
import sys

if __name__ == "__main__":
    vcf_file, chrom, pos = sys.argv[1:]
    var = vcf.Reader(filename=vcf_file)
    v = var.fetch(chrom, int(pos), int(pos)).next()
    for i in v.samples:
        print "{0}: {1}".format(i.sample, i.gt_bases)
    