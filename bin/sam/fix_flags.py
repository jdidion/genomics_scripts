#!/bin/bash
# Fix 0x02 flag for orphan reads. Requires file to be name-sorted.

import pysam

if __name__ == "__main__":
    infile = pysam.AlignmentFile(sys.argv[1], "rb")
    outfile = pysam.AlignmentFile(sys.argv[2], "wb", template=infile)
    
    try:
        l = r = None
        while True:
            if l is None:
                l = infile.next()

            r = infile.next()
        
            if l.query_name == r.query_name:
                # TODO: set 1,2
                outfile.write(l)
                outfile.write(r)
                l = r = None
                
            elif l.flag & 10 > 0:
                l.flag = l.flag & 1014 # clear bits 0x02 and 0x08
                outfile.write(l)
                l = r
                r = None
                
    except StopIteration:
        pass
    
    infile.close()
    outfile.close()


