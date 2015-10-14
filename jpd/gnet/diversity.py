import numpy as np
from math import log, exp
from itertools import product

twobit = dict((("A", 0b00), ("C", 0b01), ("G", 0b10), ("T", 0b11)))

class SeqDiv(object):
    """
    Python implementation of the DirectMap class from the SeqDiverse package written
    by Justin Brown, http://code.google.com/p/seqdiverse/. This class is used to determine
    diversity in a collection of sequences (such as HTS reads). The performance of this
    class is always linear with the number of reads (N), read length (L), and K: 
    O(N * (L - K + 1)). However, the memory requirement of this class increases exponentially
    with K. Therefore, if the read length is known ahead of time, estimate_memory() can be used
    to check that enough memory is available for the desired k.
    """
    
    @staticmethod
    def estimate_memory(read_len, k, dtype=np.int32):
        return (read_len - k + 1) * pow(4, k) * dtype().itemsize
    
    def __init__(self, k, dtype=np.int32):
        """
        The kmer counts are stored in a numpy array. The default data type of int32 should work
        for most applications but may need to be increased to int64 if there will be 
        tens-of-billions of reads or more.
        
        TODO: this implementation is way too slow for a large number of reads (can only do ~100k
        reads per minute). The core should be written in C and use Cython to populate the numpy
        array.
        """
        self.k = k
        self.nperm = pow(4, k)
        self.npos = 1
        self.data = np.zeros((self.npos, self.nperm), dtype)
    
    def _check_size(self, n):
        npos = n - self.k + 1
        if self.npos < npos:
            self.npos = npos
            self.data.resize((self.npos, self.nperm))
    
    def get_kmers(self):
        return list("".join(p) for p in product("ACGT", repeat=self.k))

    def insert(self, read):
        read_len = len(read)
        self._check_size(read_len)
        read = read.upper()
        self.insert_nocheck(read, read_len)
    
    def insert_nocheck(self, read, read_len):
        """
        This method saves a little time by not checking the read size or converting
        to upper case. Only use this method if you know that your reads are of uniform
        size and are already upper case.
        """
        
        key = 0
        check = 0
        mask = self.nperm - 1
        
        for i in xrange(0, read_len):
            b = twobit.get(read[i], 0b100)
            check = (check << 2) & mask
            if b == 0b100:
                check |= 0b11
            
            if (i >= self.k - 1) and check == 0:
                key = ((key << 2) | (b & 0b11)) & mask
                self.data[i - self.k + 1, key] += 1
    
    def write_counts(self, out):
        """Write out a matrix of counts (position x kmer) to a file object."""
        out.write(",".join(self.get_kmers()))
        out.write("\n")
        np.savetxt(out, self.data, "%i", ",")
    
    def get_summary(self):
        """
        Calculate diversity, richness and evenness at each position. Returns a numpy array
        with one row of summary statistics per position.
        """
        
        nrow = self.data.shape[0]
        summary = np.zeros((nrow, 3))
        
        for i in xrange(0, nrow):
            row = self.data[i,]
            total = sum(row)
            nonzero = row > 0
            unique = sum(nonzero)
            xlogx = sum(row[nonzero] * np.log(row[nonzero]))
            
            if total == 0:
                # consider a completely ambiguous sequence as a single species
                total = 1
                unique = 1
            
            diversity = log(total) - (xlogx / total)
            
            exp_diversity = exp(diversity) / self.nperm
            richness = float(unique) / self.nperm
            evenness = exp(diversity) / unique
            
            summary[i,] = (exp_diversity, richness, evenness)
        
        return summary
    
    def write_summary(self, out):
        """Write out a summary matrix (position x stat) to a file object."""
        summary = self.get_summary()
        out.write("diversity,richness,evenness\n")
        np.savetxt(out, summary, "%0.6f", ",")