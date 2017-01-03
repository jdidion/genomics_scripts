#!/usr/bin/env python
# Extracts read names from a SAMBLASTER discordants.bam file for which both
# reads align to chrM and either 1) read 1 is on the reverse strand and its
# position is < read 2 position, or 2) read 1 i son the forward strand and
# its position is > read 2. The reads must be within `--max-fragsize`.

from argparse import ArgumentParser
import pysam
import gzip

def main():
    parser = ArgumentParser()
    parser.add_argument("-m", "--max-fragsize", default=1000)
    parser.add_argument("-d", "--discordants", help="Must be name-sorted")
    parser.add_argument("-f", "--fastqs", nargs=2)
    parser.add_argument("-o", "--output")
    parser.add_argument("-s", "--chrM-size", type=int, default=16571)
    args = parser.parse_args()
    
    infile = pysam.AlignmentFile(args.discordants, "rb")
    outfile = pysam.AlignmentFile(args.output, "wb", template=infile)
    keepers = {}
    
    # Extract discordant reads that span the chrM artificial boundary
    for r1, r2 in zip(infile, infile):
        if r1.is_read2:
            r1, r2 = (r2, r1)
        keep = False
        if r1.reference_name == 'chrM' and r2.reference_name == 'chrM':
            if r1.is_reverse and r1.reference_start < r2.reference_start:
                keep = ((args.chrM_size - r2.reference_start) + r1.infer_query_length()) <= args.max_fragsize
            elif r2.is_reverse and r2.reference_start < r1.reference_start:
                keep = ((args.chrM_size - r1.reference_start) + r2.infer_query_length()) <= args.max_fragsize
        
        if keep:
            r1.flag = r1.flag & 16
            r1.next_reference_id = r1.next_reference_start = -1
            r1.template_length = 0
            r2.flag = r2.flag & 16
            r2.next_reference_id = r2.next_reference_start = -1
            r2.template_length = 0
            keepers[r1.query_name] = (r1,r2)
    
    # Now extract the sequence and quality from the fastq files,
    # add to the reads, and write to the output bam
    def openfq(fq):
        if fq.endswith(".gz"):
            fh = gzip.open(fq, "rt")
        else:
            fh = open(fq, "rU")
        itr = zip(*[fh] * 4)
        return (fh, itr)
    
    fq1_fh, fq1 = openfq(args.fastqs[0])
    fq2_fh, fq2 = openfq(args.fastqs[1])
    
    try:
        for r1_fq, r2_fq in zip(fq1, fq2):
            qname = r1_fq[0][1:].rstrip()
            if qname in keepers:
                r1, r2 = keepers.pop(qname)
                r1.query_sequence = r1_fq[1].rstrip().encode("utf-8")
                r1.query_qualities = "".join(chr(ord(c)-33) for c in r1_fq[3].rstrip()).encode("utf-8")
                r2.query_sequence = r2_fq[1].rstrip().encode("utf-8")
                r2.query_qualities = "".join(chr(ord(c)-33) for c in r2_fq[3].rstrip()).encode("utf-8")
                outfile.write(r1)
                outfile.write(r2)
                if len(keepers) == 0:
                    break
    finally:
        fq1_fh.close()
        fq2_fh.close()
    
    if len(keepers) > 0:
        print("Did not find sequences for {}".format(",".join(keepers.keys())))
    

if __name__ == "__main__":
    main()
