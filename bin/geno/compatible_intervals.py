#!/usr/bin/env python
# ----------------------------------------------------------------
# Given genotypes for a chromosome (or sub-chromosomal region),
# this script identifies compatible intervals, that is intervals
# within which there is no evidence of historical recombination
# among the genotyped samples. For each interval the output is
# the boundaries of the interval, the number of unique SDPs and 
# the number of unique haplotypes.
#
# Jeremy Wang
# 6/1/2010: This package is the most up-to-date and compact, readable
# code to do a complete set of compatible interval computations.
#
# 1/10/2012: Combined the critical pieces of tools.py (the beginning 
# of time) into this - it stands on its own. Trees are not yet fully 
# supported.
#
# John Didion
# 4/30/2012: Re-wrote large chunks of the code to improve performance 
# and readability. One improvement (and potential bug fix) is that
# bitarrays of homozygous calls are always relative to the minor
# allele rather than the B6 allele. This also means that the code
# is agnostic of the genotype format - if you deviate from the
# standard A/C/G/T/H/V/D/N encoding, you just need to specify the
# --missing_chars parameter. This script is no longer stand-alone.
# ----------------------------------------------------------------

from bitarray import bitarray
from collections import defaultdict, Set, OrderedDict
from csv import reader, writer
import logging
import os
import sys
import unittest

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
from util.cl import parse
from util.collections import index_map, Range

# ----------------------------------------------------------------
# Data structures
# ----------------------------------------------------------------

class SNP(object):
    def __init__(self, position, alleles, sdp_index):
        self.position = position
        self.sdp_index = sdp_index
        self.alleles = alleles
        
    def __repr__(self):
        return str(self.position)
    
    def __cmp__(self, other):
        return cmp(self.position, other.position)

class SDP:
    def __init__(self, index, hom, missing, mask):
        '''
        Each bitarray whose bits represent presence of the particular call. Each array stores the 
        call for each strain on this SDP in the order they appear in the Chromosome.strains array.
        
        hom: bitarray of homozygous calls.
        missing: dict in which keys are non-homozygous calls and values are bitarrays.
        '''
        self.index = index
        self.hom = hom
        self.missing = missing
        self.mask = mask
        
    def __len__(self):
        return len(self.hom)

    def toarray(self):
        if not hasattr(self, "_array"):
            array = list(self.hom.to01())
            for i in xrange(len(self)):
                for k,v in self.missing.iteritems():
                    if v[i]: 
                        array[i] = k
                        break
            self._array = tuple(array)
        return self._array

    def __repr__(self):
        if not hasattr(self, "_repr"):
            self._repr = ''.join(self.toarray())
        return self._repr
    
    def __eq__(self, other):
        return repr(self) == repr(other)
            
    def deepcopy(self):
        new_missing = dict((k,bitarray(v)) for k,v in self.missing.iteritems())
        return SDP(self.index, bitarray(self.hom), new_missing, bitarray(self.mask))
    
    def set_strain(self, idx, call):
        if isinstance(call, int):
            self.hom[idx] = call
            self.mask[idx] = True
            for b in self.missing.values():
                b[idx] = False
        else:
            self.mask[idx] = False
            if call not in self.missing():
                self.missing[call][idx] = True
            for k,b in self.missing.iteritems():
                b[idx] = (k == call)

def unify_haplotypes(h1, h2, method="l", collapse_chars=("N",)):
    """
    If haploytypes only differ by Ns, compress to a single haplotype.
    """
    copied = False
    for i in xrange(len(h1)):
        a = h1[i]
        b = h2[i]
        if a == b:
            continue

        if a in collapse_chars or b in collapse_chars:
            if not copied:
                h1 = list(h1)
                copied = True
            if method == 'a' and (b not in collapse_chars):
                h1[i] = b
            elif method == 'l' and (a not in collapse_chars):
                h1[i] = collapse_chars[0]
        else:
            return None

    return tuple(h1)

class Chromosome:
    def __init__(self, strains, snps, sdps):
        self.strains = strains
        self.snps = snps
        self.sdps = sdps

    def num_snps(self):
        return len(self.snps)

    def snp_positions(self, indices):
        """Returns a tuple with the genomic positions of the SNPs at the specified indices."""
        return (self.snps[i].position for i in indices)
    
    def get_all_sdps(self, start_idx=0, end_idx=None):
        """
        Returns a list of SDPs for the SNPs in the given range (inclusive). If two SNPs have
        the same SDP, then that SDP will appear twice.
        """
        if end_idx is None:
            end_idx = len(self.snps)
        return (self.sdps[s.sdp_index] for s in self.snps[start_idx:(end_idx+1)])
    
    def _get_unique_sdp_indices(self, start_idx=0, end_idx=None):
        if end_idx is None:
            end_idx = len(self.snps)
        return set((s.sdp_index for s in self.snps[start_idx:(end_idx+1)]))
        
    def count_unique_sdps(self, start_idx=0, end_idx=None):
        """Returns the number of unique SDPs in the given range (inclusive)."""
        return len(self._get_unique_sdp_indices(start_idx, end_idx))
    
    def get_unique_sdps(self, start_idx=0, end_idx=None):
        """Returns unique SDPs in the given range (inclusive)."""
        return (self.sdps[i] for i in sorted(self._get_unique_sdp_indices(start_idx, end_idx)))
    
    def get_unique_haplotypes(self, start_idx=0, end_idx=None, method="c", collapse_chars=("N",)):
        """
        Returns a map of unique haplotypes to strain indices. When there are two haplotypes that 
        differ only at sites where one haplotype is N, we consider them equal. However, this type 
        of equality is not transitive: there may be two haplotypes that are equal to a third 
        haplotype but not to each other. There are two possible ways of resolving this: 1) convert 
        Ns in each haplotype to the non-N call in the other haplotype (more conservative, 
        potentially overestimates the number of unique haplotypes), or 2) at sites where one 
        haplotype has Ns, covert the other haplotype to Ns (more liberal, potentially underestimates 
        the number of unique haplotypes).
        """
        # flatten each SDP into an array
        sdps = (s.toarray() for s in self.get_all_sdps(start_idx, end_idx)) 

        # transpose SDPs to create haplotypes and map haplotypes to strains
        hapmap = {}
        for i,h in enumerate(zip(*sdps)):
            h = tuple(h)
            if h in hapmap:
                hapmap[h].append(i)
            else:
                hapmap[h] = [i]

        # unify any haplotypes differ only by Ns
        if method != "c":
            haps = hapmap.keys()
            i = 0
            while i < len(haps) - 1:
                h1 = haps[i]
                j = i+1
                while j < len(haps):
                    h2 = haps[j]
                    new_h = unify_haplotypes(h1, h2, method, collapse_chars)
                    if new_h is None:
                        j += 1
                    elif new_h not in hapmap:
                        hapmap[new_h] = hapmap.pop(h1) + hapmap.pop(h2)
                        haps[j] = new_h
                        break
                    elif h1 == new_h:
                        hapmap[h1] += hapmap.pop(h2)
                        del haps[j]
                    elif h2 == new_h:
                        hapmap[h2] += hapmap.pop(h1)
                        break
                    else:
                        hapmap[new_h] += hapmap.pop(h1) + hapmap.pop(h2)
                        del haps[j]
                        break
                i += 1

        return hapmap
    
    def get_intervals(self, method):
        """
        Computes compatible intervals using the specified method and returns an iterable. Valid
        methods are: left, right, uber, cores and maxk.
        """
        
        # Left-to-right scan. Resulting intervals will not overlap.
        if method in ("left", "cores", "maxk"):
            if not hasattr(self, "left"):
                self.left = all_scan(self, ltor=True, uber=False)
            if method == "left": return self.left
            
        # Right-to-left scan. Resulting intervals will not overlap.
        if method in ("right", "cores", "maxk"):
            if not hasattr(self, "right"):
                self.right = all_scan(self, ltor=False, uber=False)
            if method == "right": return self.right
        
        # Uber scan (left-to-right with back-scanning). Resulting intervals may overlap.
        if method in ("uber", "maxk"):
            if not hasattr(self, "uber"):
                self.uber = all_scan(self, ltor=True, uber=True)
            if method == "uber": return self.uber
        
        # Core intervals are the intersections of corresponding intervals from left-to-right
        # and right-to-left scans.
        if method in ("cores", "maxk"):
            if not hasattr(self, "cores"):
                if len(self.left) == 0 or len(self.left) != len(self.right):
                    raise Exception("Can't core scan because left and right interval sets are "\
                        "empty or of unequal length ({0}, {1}).".format(len(self.left), len(self.right)))
                self.cores = Intervals(0, len(self.snps)-1)
                for l,r in zip(self.left, self.right):
                    self.cores.add(l.start, r.end)
            if method == "cores": return self.cores
        
        if method == "maxk":
            if not hasattr(self, "maxk"):
                self.maxk = maxk_scan(self)
            return self.maxk
        
        raise Exception("Invalid interval method: %s" % method)

# ----------------------------------------------------------------
# LR, RL, and Uber Scan computation
# ----------------------------------------------------------------

def compatible(sdp1, sdp2):
    # find all sites where neither SDP has a missing value
    mask = sdp1.mask & sdp2.mask

    # if each alleleic combination is present at at least one site, the SDPs fail the 4-gamete test
    s1 = sdp1.hom
    s2 = sdp2.hom
    s3 = bitarray(s1)
    s3.invert()
    s4 = bitarray(s2)
    s4.invert()

    return ((s1 & s2 & mask).count() == 0 or
            (s2 & s3 & mask).count() == 0 or
            (s1 & s4 & mask).count() == 0 or
            (s3 & s4 & mask).count() == 0)
            
class CompatibleTest(object):
    """
    Perform compatibility tests between SDPs. We only have to test each pair of SDPs once, so 
    we store the results of each test for future use.
    """
    def __init__(self, snps, sdps):
        self.snps = snps
        self.sdps = sdps
        self.compat_results = {}
    
    def __call__(self, snp_idx1, snp_idx2):
        sdp_idx1 = self.snps[snp_idx1].sdp_index
        sdp_idx2 = self.snps[snp_idx2].sdp_index
        
        # identical SDPs are always compatible
        if sdp_idx1 == sdp_idx2:
            return True
        
        key = tuple(sorted((sdp_idx1, sdp_idx2)))
        if key not in self.compat_results:
            sdp1 = self.sdps[sdp_idx1]
            sdp2 = self.sdps[sdp_idx2]
            self.compat_results[key] = compatible(sdp1, sdp2)
        return self.compat_results[key]

class Intervals(object):
    """
    Maintains a set of intervals. The intervals may be created in right-to-left order, in
    which case they are automatically converted to left-to-right orientation.
    """
    def __init__(self, start=0, end=1, ltor=True):
        self.start = start
        self.end = end
        self.ltor = ltor
        self.invs = []
    
    def __len__(self):
        return len(self.invs)

    def __iter__(self):
        return iter(self.invs) if self.ltor else reversed(self.invs)

    def __getitem__(self, i):
        return self.invs[i]
        
    def add(self, start, end):
        if start < self.start or end > self.end:
            raise ValueError()
        if self.ltor:
            self.invs.append(Range(start, end, (True,True)))
        else:
            self.invs.append(Range(self.end - end, self.end - start, (True,True)))
    
def all_scan(chrm, ltor=True, uber=False):
    """
    Identify intervals in which the 4-gamete test succeeds at consecutive SNPs. The ltor flag
    controls whether the scan is left-to-right (True) or right-to-left (False). If the uber flag
    is set, there is a reverse-scan from the end of the previous interval before the forward scan.
    """
    snps = chrm.snps if ltor else list(reversed(chrm.snps))
    test = CompatibleTest(snps, chrm.sdps)
    inv_start = 0
    inv_end = 1
    limit = len(snps)
    invs = Intervals(start=0, end=limit-1, ltor=ltor)
    
    while inv_end <= limit:
        next_start = scan_forward(inv_start, inv_end, limit, test)
        invs.add(inv_start, next_start - 1)
        if next_start < limit:
            if uber:
                inv_start = scan_back(next_start, test)
            else:
                inv_start = next_start
        inv_end = next_start + 1
    
    return invs

def scan_forward(inv_start, inv_end, limit, test):
    while inv_end < limit:
        # Compatibility is not transitive, so we always have to check each SDP against 
        # all other SDPs in the interval.
        for i in xrange(inv_start, inv_end):
            if not test(i, inv_end):
                return inv_end
        inv_end += 1
    return limit

def scan_back(next_start, test):
    inv_start = next_start - 1
    while inv_start >= 0:
        for i in xrange(next_start, inv_start, -1):
            if not test(i, inv_start):
                return inv_start + 1
        inv_start -= 1
    return 0

# ----------------------------------------------------------------
# MaxK data structures and computation
# ----------------------------------------------------------------

class IntervalNode:
    def __init__(self, inv):
        self.inv = inv
        self.next = []
        self.prev = None
        self.prev_count = 0
        self.total_overlap = 0
        
    def add_next(self, node):
        overlap = self.inv.end - node.inv.start + 1 
        if overlap >= 0:
            self.next.append(node)
            new_total = self.total_overlap + overlap
            if new_total == node.total_overlap:
                node.prev = self
                self.prev_count += 1
            elif new_total > node.total_overlap:
                node.total_overlap = new_total
                node.prev = self
                self.prev_count = 1
        
class IntervalLinkedList:
    def __init__(self):
        self.current = []
        self.last = None
        self.start = self.current
    
    def next_level(self):
        if len(self.current) > 0:
            self.last = self.current
            self.current = []
        
    def add_node(self, inv):
        new_node = IntervalNode(inv)
        self.current.append(new_node)
        if self.last is not None:
            for node in self.last:
                node.add_next(new_node)
                
    def last_node(self):
        max_total = -1
        last = None
        for node in self.current:
            if node.total_overlap >= max_total:
                max_total = node.total_overlap
                last = node
        return last
    
def maxk_scan(chrm):
    uber = chrm.get_intervals("uber")
    core = chrm.get_intervals("cores")

    min_list = IntervalLinkedList()
    last = 0
    for i in xrange(len(core)):
        min_list.next_level()
        for j in xrange(last, len(uber)):
            if uber[j].start > core[i].start:
                break
            elif (uber[j].end >= core[i].end):
                min_list.add_node(uber[j])
        last = j
        
    invs = []
    node = min_list.last_node()
    while node != None:
        invs.append(node.inv)
        node = node.prev
    invs.reverse()
    return invs

# ----------------------------------------------------------
# Tree computation
# ----------------------------------------------------------

def simple_dist(h1, h2):
    """
    A simple method of distance calculation: [number of sites that are homozygous in both
    haplotypes and unequal] / [total number of sites that are homozygous in both haplotypes].
    """
    diff = 0
    total = 0
    for a,b in zip(h1,h2):
        if a in ('0','1') and b in ('0','1'):
            total += 1
            if a != b:
                diff += 1
    if total == 0:
        return 0
    else:
        return float(diff) / total

def compute_tree(haps, dist_fn=simple_dist):
    n = len(haps)
    
    # no need to compute a tree if we only have one or two haplotypes
    if n == 1:
        return "(0);"
    elif n == 2:
        return "(0,1):{0}".format(simple_dist(haps[0], haps[1]))
    
    dists = {}
    
    for i in xrange(n):
        for j in xrange(i+1, n):
            dists[(i,j)] = dist_fn(haps[i], haps[j])
    
    (result,) = nj.gnj(dists, keep=1, show_progress=False)
    (score, tree) = result
    return str(tree)

# ----------------------------------------------------------
# I/O
# ----------------------------------------------------------

class InfileReader(object):
    """
    Read an input file where the first row is a header and two column indices determine
    which columns will be used: genomic position and first genotype column. The reader stores 
    the strain list (first row, starting from the first genotype column) and iterates over 
    subsequent rows. Each row will be returned as a two-element tuple:
    ( pos, (geno1, geno2, ...) ).
    """
    def __init__(self, path, delim, columns):
        self.handle = open(path, 'rU')
        self.reader = reader(self.handle, delimiter=delim)
        self.strains = index_map(self.reader.next()[columns[1]:])
        self._pos_column = columns[0]
        self._first_geno_column = columns[1]
        
    def __iter__(self):
        return self

    def next(self):
        try:
            row = self.reader.next()
            return (int(row[self._pos_column]), row[self._first_geno_column:])
        except StopIteration:
            self.handle.close()
            raise

def parse_input(geno_reader, max_n_frac=1.0, ignore_file=None, missing_chars=('N','H','V','D')):
    N_char = missing_chars[0]
    missing_chars = set(missing_chars)
    length = len(geno_reader.strains)
    max_ns = length * max_n_frac
    
    sdp_map = OrderedDict()
    num_sdps = 0
    snps = []
    num_lines = 0
    not_informative = 0
    too_many_alleles = 0
    too_many_ns = 0
    
    init = '0' * length
    def make_bitarray():
        return bitarray(init)
    
    for s in geno_reader:
        num_lines += 1
        pos = s[0]
        
        # Determine the A and B alleles. Maintain a bitarray for homozygous alleles where 0=A
        # and 1=B. Also maintain a separate bitarray for each N character.
        homozygous = defaultdict(make_bitarray) # homozygous samples
        missing = defaultdict(make_bitarray) # heterozygous or other samples
        
        for i,g in enumerate(s[1]):
            if g in missing_chars:
                missing[g][i] = True
            else:
                homozygous[g][i] = True
        
        if len(missing) > 0:
            ns = reduce(lambda x,y: x|y, missing.values())

            # skip this SNP if there are too many Ns
            if ns.count() > max_ns:
                too_many_ns += 1
                continue
        
            # mask is a bitarray representing the indices of SNPs with homozygous calls
            mask = ~ns
        else:
            mask = bitarray('1' * length)

        # order alleles by increasing frequency
        freq = dict((k, v.count()) for k,v in homozygous.iteritems())
        alleles = sorted(freq.keys(), key=lambda k: freq[k])
        
        if len(alleles) < 2:
            # represent a non-informative bitarray as all zeros
            not_informative += 1
            hom = make_bitarray()
        else:
            # if there are more than two alleles, convert the most minor to Ns
            if len(alleles) > 2:
                too_many_alleles += 1
                num_remove = len(alleles) - 2
                for a in alleles[0:num_remove]:
                    missing[N_char] |= homozygous[a]
                alleles = alleles[num_remove:]
                
            # if two alleles have equal frequency, choose the one with the first set bit
            if freq[alleles[0]] < freq[alleles[1]]:
                hom = homozygous[alleles[0]]
            else:
                for i in xrange(0, length):
                    if homozygous[alleles[0]][i]:
                        hom = homozygous[alleles[0]]
                        break
                    elif homozygous[alleles[1]][i]:
                        hom = homozygous[alleles[1]]
                        alleles = reversed(alleles)
                        break
        
        sdp = SDP(num_sdps, hom, missing, mask)
        sdp_key = str(sdp)
        
        if sdp_key not in sdp_map:
            sdp_map[sdp_key] = sdp
            num_sdps += 1
        else:
            sdp = sdp_map[sdp_key]
        
        snp = SNP(pos, alleles, sdp.index)
        snps.append(snp)
    
    logging.info("Total lines: {0}".format(num_lines))
    logging.info("Too many N's: {0}".format(too_many_ns))
    logging.info("Total SNPs: {0}".format(len(snps)))
    logging.info("<2 alleles: {0}".format(not_informative))
    logging.info(">2 alleles: {0}".format(too_many_alleles))
    
    if ignore_file != None:
        with open(ignore_file, 'rU') as ifile:
            for row in reader(ifile):
                if row[0] not in geno_reader.strains:
                    continue
                
                strain_idx = geno_reader.strains[row[0]]
                start = int(row[1])
                end = int(row[2])
                call = row[3]
            
                while s < len(snps):
                    if snps[s].position > end:
                        break
                    elif snps[s].position >= start:
                        # Make a copy of the SDP, convert the specified strain to the specified
                        # call and update the SNP reference. We don't modify the SDP in place
                        # since other SNPs may reference it.
                        sdp = sdps[snps[s].sdp_index].deepcopy()
                        sdp.index = num_sdps
                        sdp.set_strain(strain_idx, call)
                        sdp_key = repr(sdp)
                        if sdp_key not in sdp_map:
                            sdp_map[sdp_key] = sdp
                            num_sdps += 1
                        else:
                            sdp = sdp_map[sdp_key]
                        snps[s].sdp_index = sdp.index
    
    return Chromosome(geno_reader.strains, snps, sdp_map.values())

def write_output(chrm, invs, outfile, remove_singletons=False, haplotype_method="c", collapse_chars=("N",),
        write_haplotypes=False, write_trees=False):
    """Writes a tab-delimited file in which each row is an interval."""
    
    with open(outfile, 'w') as o:
        w = writer(o, delimiter="\t")
        
        header = ["minstart", "minend","midstart", "midend", "maxstart",  "maxend", 
                  "minsize", "midsize", "maxsize", "nsnps", "nsdps", "nhaps"]
        if write_haplotypes:
            header.append("haps")
        if write_trees:
            header.extend(("tree", "samples"))
        w.writerow(header)
        
        num_singletons = 0
        
        for i,inv in enumerate(invs):
            snp_count = len(inv)
            if remove_singletons and snp_count == 1:
                num_singletons += 1
                continue
                
            start = inv.start
            end = inv.end
            
            # The minimum range is between the SNP after the last SNP of the previous interval
            # and the SNP before the first SNP of the next interval
            min_rng = list(chrm.snp_positions((
                start if i == 0 else invs[i-1].end + 1,
                end if i == len(invs)-1 else invs[i+1].start - 1)))
            # The middle range is between the first and last SNP in the interval
            mid_rng = list(chrm.snp_positions((start, end)))
            # The outer range is between the end of the previous interval and the beginning of
            # the next interval.
            max_adj = (0 if start == 0 else 1, 0 if end == chrm.num_snps() - 1 else 1)
            max_rng = list(chrm.snp_positions((start - max_adj[0], end + max_adj[1])))
            max_rng[0] += max_adj[0]
            max_rng[1] -= max_adj[1]
            
            haps = chrm.get_unique_haplotypes(start, end, haplotype_method, collapse_chars)
            
            row = min_rng + mid_rng + max_rng
            row.append(min_rng[1] - min_rng[0] + 1) # inner_size
            row.append(mid_rng[1] - mid_rng[0] + 1) # middle_size
            row.append(max_rng[1] - max_rng[0] + 1) # outer_size
            row.append(snp_count)
            row.append(chrm.count_unique_sdps(start, end))
            row.append(len(haps))
            
            if write_haplotypes:
                row.append('|'.join(('='.join((str(i), ''.join(h))) for i,h in enumerate(haps.keys()))))
                
            if write_trees:
                # append the newick format of the haplotype tree and the set of
                # strains associated with each branch
                row.append(compute_tree(haps.keys()))
                row.append('|'.join(','.join(str(s) for s in g) for g in haps.values()))
            
            w.writerow(row)
        
    if remove_singletons:
        logging.info("Number of singletons removed: {0}".format(num_singletons))

if __name__ == "__main__":
    def add_args(parser):
        parser.add_argument("-c", "--columns", metavar="INDEX", type='int_list', default=[1,2],
            help="List of two column indicies (0-based): position and first genotype.")
        parser.add_argument("-d", "--delim", metavar="CHAR", default=",",
            help="Field delimiter in input and output files.")
        parser.add_argument("-H", "--haplotype_method", metavar="METHOD", choices=["c","a","l"], default="a",
            help="Method for identifying unique haplotypes. c=only remove exact duplicates (most "\
                 "conservative), a=convert Ns to non-N calls (arbitrary), l=convert non-N calls to "\
                 "Ns (most liberal)")
        parser.add_argument("-i", "--ignore_file", type="readable_file", metavar="FILE", default=None,
            help="File containing regions to ignore.")
        parser.add_argument("-m", "--missing_chars", type="str_list", metavar="LIST", default=("N","H","V","D"),
            help="List of characters to treat as missing information. The first character should "\
                 "be the 'no-call' character.")
        parser.add_argument("-n", "--max_n", metavar="FRACTION", type=float, default=0.2,
            help="Maximum fraction of N calls for a SNP.")
        parser.add_argument("-o", "--collapse_missing", action="store_true", default=False,
            help="Whether to collapse haplotypes that differ only by missing characters. The default is to "\
                 "collapse only haplotypes that differ by N-calls.")
        parser.add_argument("-t", "--interval_type", metavar="TYPE", 
            choices=["left", "right", "cores", "uber", "maxk"], default="maxk", 
            help="Algorithm for determining intervals.")
        parser.add_argument("--remove_singletons", action="store_true", default=False,
            help="Remove single-SNP intervals from output.")
        parser.add_argument("--write_haplotypes", action="store_true", default=False,
            help="Write unique haplotypes for each interval to the output file.")
        parser.add_argument("--write_trees", action="store_true", default=False,
            help="Compute and output a distance tree for each interval.")
        parser.add_argument("infile", type="readable_file", metavar="FILE",
            help="File contaning genotypes.")
        parser.add_argument("outfile", type="writeable_file", metavar="FILE",
            help="Output file.")
            
    ns = parse(add_args)
    
    logging.debug("Parsing input file")
    reader = InfileReader(ns.infile, ns.delim, ns.columns)
    chrm = parse_input(reader, ns.max_n, ns.ignore_file, ns.missing_chars)
    
    logging.debug("Computing intervals")
    invs = chrm.get_intervals(ns.interval_type)
    logging.info("Number of intervals: {0}".format(len(invs)))
    
    logging.debug("Writing output")
    if ns.write_trees:
        # This is a hack to avoid warnings from PyCogent about not having MPI enabled
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from cogent.phylo import nj
    
    if args.collapse_missing:
        collapse_chars = ns.missing_chars
    else:
        collapse_chars = (ns.missing_chars[0],)
    write_output(chrm, invs, ns.outfile, ns.remove_singletons, ns.haplotype_method, 
        collapse_chars, ns.write_haplotypes, ns.write_trees)

class MockReader(object):
    def __init__(self, rows, header=None):
        self.rows = rows
        if header is None:
            header = (str(x) for x in xrange(0, len(rows[0][1])))
        self.strains = index_map(header)
        
    def __iter__(self):
        return iter(self.rows)

class CompatInvTests(unittest.TestCase):
    def test_parse_simple(self):
        reader = MockReader([(1, ('A','C','A','N'))])
        chrm = parse_input(reader)
        self.assertIsNotNone(chrm)
        self.assertEquals(chrm.num_snps(), 1)
        self.assertEquals(chrm.count_unique_sdps(), 1)
        self.assertEquals(chrm.snps[0].position, 1)
        self.assertEquals(chrm.snps[0].sdp_index, 0)
        self.assertEquals(chrm.snps[0].alleles, ['C','A'])
        sdps = list(chrm.get_all_sdps())
        self.assertEquals(len(sdps), 1)
        self.assertEquals(sdps[0].index, 0)
        self.assertEquals(sdps[0].hom, bitarray('0100'))
        self.assertEquals(len(sdps[0].missing), 1)
        self.assertEquals(sdps[0].missing['N'], bitarray('0001'))
        self.assertEquals(sdps[0].mask, bitarray('1110'))
        self.assertEquals(str(sdps[0]), '010N')

    def test_parse_equal_maf(self):
        reader = MockReader([(1, ('A','C','A','C'))])
        chrm = parse_input(reader)
        self.assertEquals(chrm.snps[0].alleles, ['A','C'])
        self.assertEquals(chrm.sdps[0].hom, bitarray('1010'))
        self.assertEquals(chrm.sdps[0].mask, bitarray('1111'))
    
    def test_parse_too_many_ns(self):
        reader = MockReader([(1, ('A','N','A','N'))])
        chrm = parse_input(reader, max_n_frac=0.2)
        self.assertEquals(chrm.num_snps(), 0)
        self.assertEquals(chrm.count_unique_sdps(), 0)
        
    def test_three_alleles(self):
        reader = MockReader([(1, ('A','C','A','T','C','A'))])
        chrm = parse_input(reader)
        self.assertEquals(chrm.snps[0].alleles, ['C','A'])
        self.assertEquals(str(chrm.sdps[0]), '010N10')
        
    def test_same_sdp(self):
        reader = MockReader([
            (1, ('A','N','A','N','C')),
            (2, ('C','N','C','N','A'))
        ])
        chrm = parse_input(reader)
        self.assertEquals(chrm.num_snps(), 2)
        self.assertEquals(chrm.count_unique_sdps(), 1)
        self.assertEquals(chrm.sdps[0].hom, bitarray('00001'))
        self.assertEquals(chrm.snps[0].sdp_index, 0)
        self.assertEquals(chrm.snps[1].sdp_index, 0)
        
    def test_same_compatible(self):
        reader = MockReader([
            (1, ('C','A')),
            (2, ('C','A'))
        ])
        chrm = parse_input(reader)
        # identical SDPs should be compatible
        sdps = list(chrm.get_all_sdps())
        self.assertEquals(len(sdps), 2)
        self.assertEquals(sdps[0].index, sdps[1].index)
        self.assertTrue(compatible(sdps[0], sdps[1]))
        
        reader = MockReader([
            (1, ('C','A','C')),
            (2, ('C','A','N'))
        ])
        chrm = parse_input(reader)
        # SDPs that differ only by missing values should also be compatible
        sdps = list(chrm.get_all_sdps())
        self.assertEquals(len(sdps), 2)
        self.assertTrue(compatible(sdps[0], sdps[1]))

    def test_not_compatible(self):
        reader = MockReader([
            (1, ('C','A','C','A')),
            (2, ('A','C','C','A'))
        ])
        chrm = parse_input(reader)
        sdps = list(chrm.get_all_sdps())
        self.assertEquals(len(sdps), 2)
        self.assertNotEquals(sdps[0].index, sdps[1].index)
        self.assertFalse(compatible(sdps[0], sdps[1]))
    
    def test_intervals_ltor(self):
        i = Intervals(0, 10)
        i.add(0, 5)
        i.add(6, 10)
        invs = list(iter(i))
        self.assertEquals(range(0,6), list(iter(invs[0])))
        self.assertEquals(range(6,11), list(iter(invs[1])))
    
    def test_intervals_rtol(self):
        i = Intervals(0, 10, ltor=False)
        i.add(0, 5)
        i.add(6, 10)
        invs = list(iter(i))
        self.assertEquals(range(0,5), list(iter(invs[0])))
        self.assertEquals(range(5,11), list(iter(invs[1])))
    
    def test_get_intervals_left(self):
        reader = MockReader([
            (1, ('C','A','C','A')),
            (2, ('C','C','C','C')),
            (3, ('C','A','A','C')),
            (4, ('A','A','A','A'))
        ])
        chrm = parse_input(reader)
        self.assertEquals(4, chrm.num_snps())
        ltor = list(iter(chrm.get_intervals("left")))
        self.assertEquals(2, len(ltor))
        self.assertEquals(range(0,2), list(iter(ltor[0])))
        self.assertEquals(range(2,4), list(iter(ltor[1])))

    def test_get_intervals_right(self):
        reader = MockReader([
            (1, ('C','A','C','A')),
            (2, ('C','C','C','C')),
            (3, ('C','A','A','C')),
            (4, ('A','A','A','A'))
        ])
        chrm = parse_input(reader)
        self.assertEquals(4, chrm.num_snps())
        ltor = list(iter(chrm.get_intervals("right")))
        self.assertEquals(2, len(ltor))
        self.assertEquals(range(0,1), list(iter(ltor[0])))
        self.assertEquals(range(1,4), list(iter(ltor[1])))
        
    def test_get_intervals_uber(self):
        reader = MockReader([
            (1, ('C','A','C','A')),
            (2, ('C','C','C','C')),
            (3, ('C','A','A','C')),
            (4, ('A','A','A','A'))
        ])
        chrm = parse_input(reader)
        self.assertEquals(4, chrm.num_snps())
        ltor = list(iter(chrm.get_intervals("uber")))
        self.assertEquals(2, len(ltor))
        self.assertEquals(range(0,2), list(iter(ltor[0])))
        self.assertEquals(range(1,4), list(iter(ltor[1])))

    def test_get_intervals_cores(self):
        reader = MockReader([
            (1, ('C','A','C','A')),
            (2, ('C','C','C','C')),
            (3, ('C','A','A','C')),
            (4, ('A','A','A','A'))
        ])
        chrm = parse_input(reader)
        self.assertEquals(4, chrm.num_snps())
        ltor = list(iter(chrm.get_intervals("cores")))
        self.assertEquals(2, len(ltor))
        self.assertEquals(range(0,1), list(iter(ltor[0])))
        self.assertEquals(range(2,4), list(iter(ltor[1])))

    def test_unify_haplotypes_arbitrary(self):
        hap1 = ('0','1','N')
        hap2 = ('N','1','0')
        self.assertEquals(('0','1','0'), unify_haplotypes(hap1, hap2, 'a'))
        
    def test_unify_haplotypes_liberal(self):
        hap1 = ('0','1','N')
        hap2 = ('N','1','0')
        self.assertEquals(('N','1','N'), unify_haplotypes(hap1, hap2, 'l'))

    def test_get_unique_haplotypes_conservative(self):
        reader = MockReader([
            (1, ('A','N','C')),
            (2, ('A','A','A'))
        ])
        chrm = parse_input(reader)
        haps = chrm.get_unique_haplotypes()
        self.assertEquals(3, len(haps))

    def test_get_unique_haplotypes_arbitrary(self):
        reader = MockReader([
            (1, ('A','N','C')),
            (2, ('A','A','A'))
        ])
        chrm = parse_input(reader)
        haps = chrm.get_unique_haplotypes(method='a')
        self.assertEquals(2, len(haps))

    def test_get_unique_haplotypes_liberal(self):
        reader = MockReader([
            (1, ('A','N','C')),
            (2, ('A','A','A'))
        ])
        chrm = parse_input(reader)
        haps = chrm.get_unique_haplotypes(method='l')
        self.assertEquals(1, len(haps))
    
    def test_simple_dist(self):
        hap1 = ('0','1','0','N')
        hap2 = ('1','N','0','1')
        self.assertEquals(0.5, simple_dist(hap1, hap2))