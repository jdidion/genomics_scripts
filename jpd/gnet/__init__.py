from __future__ import absolute_import
from collections import namedtuple, OrderedDict

from jpd.util.collections import BitfieldEnum
from jpd.util.io import safe_read_file_array

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

def karyotype(path, bands=False):
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
    return (tuple(alleles), count_alleles)
