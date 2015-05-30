#!/usr/bin/env python
# Fill in missing genotypes with the reference allele following a merge.

from vcf import Reader, Writer
import vcf.model
import sys

def cumsum(it):
    total = 0
    for x in it:
        total += x
        yield total

def main():
    close_infile = False
    infile = sys.argv[1]
    if infile == "-":
        i = sys.stdin
    else:
        i = open(infile, "rU")
        close_infile = True
        
    close_outfile = False
    outfile = sys.argv[2]
    if outfile == "-":
        o = sys.stdout
    else:
        o = open(outfile, "w")
        close_outfile = True
        
    group_sizes = tuple(int(i) for i in sys.argv[3].split(","))
    group_cnt = range(len(group_sizes))
    group_ixs = (0,) + tuple(cumsum(group_sizes))
    
    reader = Reader(i)
    writer = Writer(o, reader)
    for rec in reader:
        fix = []

        for i in group_cnt:
            calls = rec.samples[group_ixs[i]:group_ixs[i+1]]
            called = False
            for c in calls:
                if c.called:
                    called = True
                    break
            if not called:
                fix.extend(calls)
        
        if len(fix) > 0:
            for c in fix:
                # This is a hack because PyVCF _Call objects are not mutable
                c.data = vcf.model.make_calldata_tuple(c.data._fields)(GT="0/0", DP=c.data.DP, GQ=c.data.GQ, PL=c.data.PL)
        
        writer.write_record(rec)
    
    if close_infile:
        i.close()
    if cloes_outfile:
        o.close()

if __name__ == "__main__":
    main()