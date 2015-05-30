from __future__ import absolute_import
from collections import namedtuple, OrderedDict

from .collections import BitfieldEnum
from .io import safe_read_file_array

# nucleotide and ambiguity codes
nt = BitfieldEnum('A','C','G','T')
nt.add_or('R', 'A', 'G') # purine
nt.add_or('Y', 'C', 'T') # pyrimidine
nt.add_or('K', 'G', 'T') # keto
nt.add_or('M', 'A', 'C') # amino
nt.add_or('S', 'C', 'G') # strong
nt.add_or('W', 'A', 'T') # weak
nt.add_or('B', 'C', 'G', 'T') # not A
nt.add_or('D', 'A', 'G', 'T') # not C
nt.add_or('H', 'A', 'C', 'T') # not G
nt.add_or('V', 'A', 'C', 'G') # not T (U is Uracil)
nt.add_or('N', 'A', 'C', 'G', 'T')

complement = dict(A='T', C='G', G='C', T='A')
def revcomp(seq):
    return ''.join(complement[c] for c in reversed(seq))

Chromosome = namedtuple('Chromosome', ['name','size','display','id','index'])
Band = namedtuple('Band', ['name','start','end','info'])

def sort_chrm(x, y):
    if isinstance(x, Chromosome): x = x.name
    if x.isdigit(): x = int(x)
    if isinstance(y, Chromosome): y = y.name
    if y.isdigit(): y = int(y)
    return cmp(x, y)

def chromosomes(*args):
    return OrderedDict((k.name, k) for k in karyotype(*args))

def karyotype(path=None, bands=False):
    if path is None:
        path = '/Volumes/Data/data/genomes/karyotype.mouse.mm9.txt'
    chrms = []
    bands = {} if bands else None
    chrm_idx = 1
    
    for line in safe_read_file_array(path):
        fields = line.split(" ")
        if fields[0] == 'chr':
            chrms.append(Chromosome(fields[3], int(fields[5]), fields[6], fields[2], chrm_idx))
            chrm_idx += 1
        elif bands and fields[0] == 'band':
            if fields[1] not in bands:
                bands[fields[1]] = []
            bands[fields[1]].append(*fields[2:6])
        else:
            break
    
    return (chrms,bands) if bands else chrms

def count_alleles(chars, N_chars=set('N')):
    alleles = set()
    counts = {}
    for a in chars:
        if a not in N_chars:
            alleles.add(a)
        if a not in counts:
            counts[a] = 0
        counts[a] += 1
    return (tuple(alleles), counts)

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