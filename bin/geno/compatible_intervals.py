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
#
# 7/12/2015
# Added ability to read genotypes from PLINK binary files.
# ----------------------------------------------------------------

from bitarray import bitarray
from collections import defaultdict, Set, OrderedDict
from csv import reader, writer
import logging
import os
import sys
import unittest

from jpd.util.cl import parse
from jpd.util.collections import index_map, Range

# ----------------------------------------------------------------
# Data structures
# ----------------------------------------------------------------

class SNP(object):
    def __init__(self, snpID, position, alleles, sdp_index):
        self.snpID = snpID
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
        call for each sample on this SDP in sample order.
        
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
            self._repr = ''.join(str(i) for i in self.toarray())
        return self._repr
    
    def __eq__(self, other):
        return repr(self) == repr(other)
            
    def deepcopy(self):
        new_missing = dict((k,bitarray(v)) for k,v in self.missing.iteritems())
        return SDP(self.index, bitarray(self.hom), new_missing, bitarray(self.mask))
    
    def set_sample(self, idx, call):
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
    def __init__(self, snps, sdps):
        self.snps = snps
        self.sdps = sdps

    def num_snps(self):
        return len(self.snps)

    def snp_ids(self, indices):
        """Returns a tuple with IDs of the SNPs at the specified indices."""
        return (self.snps[i].snpID for i in xrange(indices[0], indices[1]+1))

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
        Returns a map of unique haplotypes to sample indices. When there are two haplotypes that 
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

        # transpose SDPs to create haplotypes and map haplotypes to samples
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
    
    def get_intervals(self, method, **kwargs):
        """
        Computes compatible intervals using the specified method and returns an iterable. Valid
        methods are: left, right, uber, cores and maxk.
        """
        
        # Left-to-right scan. Resulting intervals will not overlap.
        if method in ("left", "cores", "maxk"):
            if not hasattr(self, "left"):
                self.left = all_scan(self, ltor=True, uber=False, **kwargs)
            if method == "left": return self.left
            
        # Right-to-left scan. Resulting intervals will not overlap.
        if method in ("right", "cores", "maxk"):
            if not hasattr(self, "right"):
                self.right = all_scan(self, ltor=False, uber=False, **kwargs)
            if method == "right": return self.right
        
        # Uber scan (left-to-right with back-scanning). Resulting intervals may overlap.
        if method in ("uber", "maxk"):
            if not hasattr(self, "uber"):
                self.uber = all_scan(self, ltor=True, uber=True, **kwargs)
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
    def __init__(self, snps, sdps, max_snps=None, max_size=None, max_spacing=None):
        self.snps = snps
        self.sdps = sdps
        self.max_snps = max_snps
        self.max_size = max_size
        self.max_spacing = max_spacing
        self._compat_results = {}
    
    def __call__(self, snp_idx1, snp_idx2):
        if self.max_snps is not None and abs(snp_idx2 - snp_idx1) + 1 > self.max_snps:
            return False
        
        snp_idx1, snp_idx2 = sorted((snp_idx1, snp_idx2))
        
        snp1 = self.snps[snp_idx1]
        snp2 = self.snps[snp_idx2]
        
        if self.max_size is not None and snp2.position - snp1.position + 1 > self.max_size:
            return False
        
        if self.max_spacing is not None and (
                self.snps[snp_idx1+1].position - snp1.position + 1 > self.max_spacing or (
                    snp_idx2 - snp_idx1 > 1 and 
                    snp2.position - self.snps[snp_idx2[1]].position + 1 > self.max_spacing  
                )
            ): return False
        
        sdp_idx1 = snp1.sdp_index
        sdp_idx2 = snp2.sdp_index
        
        # identical SDPs are always compatible
        if sdp_idx1 == sdp_idx2:
            return True
        
        key = tuple(sorted((sdp_idx1, sdp_idx2)))
        if key not in self._compat_results:
            sdp1 = self.sdps[sdp_idx1]
            sdp2 = self.sdps[sdp_idx2]
            self._compat_results[key] = compatible(sdp1, sdp2)

        return self._compat_results[key]

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
    
def all_scan(chrm, ltor=True, uber=False, **kwargs):
    """
    Identify intervals in which the 4-gamete test succeeds at consecutive SNPs. The ltor flag
    controls whether the scan is left-to-right (True) or right-to-left (False). If the uber flag
    is set, there is a reverse-scan from the end of the previous interval before the forward scan.
    
    max_snps: maximum size of an interval in number of SNPs
    max_size: maximum size of an interval in bp
    max_spacing: maximum distance between consecutive SNPs in bp
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
    Read an input file where the first row is a header and four column indices determine
    which columns will be used: snpID, chromosome, genomic position and first genotype column. 
    The reader stores the sample list (first row, starting from the first genotype column) 
    and iterates over subsequent rows. Each row is returned as a two-element tuple: 
    ( name, pos, (geno1, geno2, ...) ).
    """
    def __init__(self, path, chrm=None, delim=",", columns=(0,1,2,3)):
        self.chrm = chrm
        self._handle = open(path, 'rU')
        self._reader = reader(self.handle, delimiter=delim)
        self._samples = index_map(self._reader.next()[columns[1]:])
        self._id_column = columns[0]
        self._chrm_column = columns[1]
        self._pos_column = columns[2]
        self._first_geno_column = columns[3]
    
    @property
    def num_samples(self):
        return len(self._samples)
    
    def has_sample(self, name):
        return name in self._samples
    
    def get_sample_index(self, name):
        return self._samples[sample]
    
    def __iter__(self):
        return self

    def next(self):
        try:
            row = self._reader.next()
            while self.chrm is not None and str(row[self._chrm_column]) != self.chrm:
                row = self._reader.next()
            return (
                row[self._id_column], 
                int(row[self._pos_column]), 
                row[self._first_geno_column:]
            )
        except StopIteration:
            self.close()
            raise
    
    def close(self):
        self._handle.close()

class PlinkReader(object):
    """
    Iterate over data from PLINK binary files (BED/BIM/FAM). Each row is returned as a 
    two-element tuple: ( pos, (geno1, geno2, ...) ).
    """
    def __init__(self, path, chrm=None, families=None):
        from plinkio import plinkfile
        from itertools import izip
        
        self.chrm = chrm
        self.handle = plinkfile.open(path)
        if not self.handle.one_locus_per_row():
            raise Exception("This script requires that SNPs are rows and samples columns.")
        
        samples = self.handle.get_samples()
        self._subset_idxs = None
        if families is not None:
            families = set(families)
            self._subset_idxs = set(i for i,sample in enumerate(samples) if sample.fid in families)
            self._samples = index_map(samples[i].iid for i in self._subset_idxs)
        
        else:
            self._samples = dict((s.iid, s) for s in samples)
        
        self._loci = self.handle.get_loci()
        self._iter = izip(self._loci, self.handle)
    
    @property
    def num_samples(self):
        return len(self._samples)
    
    def has_sample(self, name):
        return name in self._samples
    
    def get_sample_index(self, name):
        return self._samples[sample]
        
    def __iter__(self):
        return self
    
    def next(self):
        try:
            locus, genotypes = self._iter.next()
            while self.chrm is not None and str(locus.chromosome) != self.chrm:
                locus, genotypes = self._iter.next()
            if self._subset_idxs is None:
                genotypes = list(str(g) for g in genotypes)
            else:
                genotypes = list(str(g) for i,g in enumerate(genotypes) if i in self._subset_idxs)
            return (locus.name, locus.bp_position, genotypes)
        except StopIteration:
            self.close()
            raise
    
    def close(self):
        self.handle.close()

def parse_input(geno_reader, min_maf=None, max_n_frac=1.0, ignore_file=None, 
        het_chars=('H',), missing_chars=('N','V','D')):
    N_char = missing_chars[0]
    het_chars = set(het_chars)
    missing_chars = set(missing_chars)
    length = geno_reader.num_samples
    max_ns = length * max_n_frac
    
    sdp_map = OrderedDict()
    num_sdps = 0
    snps = []
    num_lines = 0
    not_informative = 0
    too_many_alleles = 0
    too_few_hom = 0
    too_many_ns = 0
    too_low_maf = 0
    
    init = '0' * length
    def make_bitarray():
        return bitarray(init)
    
    for name, pos, genotypes in geno_reader:
        num_lines += 1
        
        # Determine the A and B alleles. Maintain a bitarray for homozygous alleles where 0=A
        # and 1=B. Also maintain a separate bitarray for each N character.
        homozygous = defaultdict(make_bitarray) # homozygous samples
        missing = defaultdict(make_bitarray) # heterozygous or other samples
        het_count = 0
        missing_count = 0
        
        for i,g in enumerate(genotypes):
            if g in het_chars:
                het_count += 1
                missing[g][i] = True
            elif g in missing_chars:
                missing_count += 1 
                missing[g][i] = True
            else:
                homozygous[g][i] = True
        
        # TODO should probably have an option for min homozygous ratio
        if len(homozygous) == 0:
            too_few_hom += 1
            continue
        
        if missing_count > max_ns:
            too_many_ns += 1
            continue
        
        if len(missing) > 0:
            # mask is a bitarray representing the indices of SNPs with homozygous calls
            ns = reduce(lambda x,y: x|y, missing.values())
            mask = ~ns
        else:
            mask = bitarray('1' * length)

        # compute allele frequency and check MAF
        freq = dict((k, v.count()) for k,v in homozygous.iteritems())
        if min_maf is not None:
            maf = float((min(freq.values()) * 2) + het_count) / float(2 * (sum(freq.values()) + het_count))
            if maf < min_maf:
                too_low_maf += 1
                continue
        
        # order alleles by increasing frequency
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
        
        snp = SNP(name, pos, alleles, sdp.index)
        snps.append(snp)
    
    logging.info("Total lines: {0}".format(num_lines))
    logging.info("Too few homozygous: {0}".format(too_few_hom))
    logging.info("Too many N's: {0}".format(too_many_ns))
    if min_maf is not None:
        logging.info("MAF < {0}: {1}".format(min_maf, too_low_maf))
    logging.info("Total SNPs: {0}".format(len(snps)))
    logging.info("<2 alleles: {0}".format(not_informative))
    logging.info(">2 alleles: {0}".format(too_many_alleles))
    
    if ignore_file != None:
        with open(ignore_file, 'rU') as ifile:
            for row in reader(ifile):
                if not geno_reader.has_sample(row[0]):
                    continue
                
                sample_idx = geno_reader.get_sample_index(row[0])
                start = int(row[1])
                end = int(row[2])
                call = row[3]
            
                for snp in snps:
                    if snp.position > end:
                        break
                    
                    elif snp.position >= start:
                        # Make a copy of the SDP, convert the specified sample to the specified
                        # call and update the SNP reference. We don't modify the SDP in place
                        # since other SNPs may reference it.
                        sdp = sdps[snp.sdp_index].deepcopy()
                        sdp.index = num_sdps
                        sdp.set_sample(sample_idx, call)
                        sdp_key = repr(sdp)
                        if sdp_key not in sdp_map:
                            sdp_map[sdp_key] = sdp
                            num_sdps += 1
                        else:
                            sdp = sdp_map[sdp_key]
                        snp.sdp_index = sdp.index
    
    return Chromosome(snps, sdp_map.values())

def write_output(chrm, invs, outfile, remove_singletons=False, haplotype_method="c", collapse_chars=("N",),
        write_haplotypes=False, write_trees=False):
    """Writes a tab-delimited file in which each row is an interval."""
    
    with open(outfile, 'w') as o:
        w = writer(o, delimiter="\t")
        
        header = ["minstart", "minend","midstart", "midend", "maxstart",  "maxend", 
                  "minsize", "midsize", "maxsize", "nsnps", "snpIDs", "nsdps", "nhaps"]
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
            snpIDs = list(chrm.snp_ids((start, end)))
            row.append(",".join(snpIDs))
            row.append(chrm.count_unique_sdps(start, end))
            row.append(len(haps))
            
            if write_haplotypes:
                row.append('|'.join(('='.join((str(i), ''.join(h))) for i,h in enumerate(haps.keys()))))
                
            if write_trees:
                # append the newick format of the haplotype tree and the set of
                # samples associated with each branch
                row.append(compute_tree(haps.keys()))
                row.append('|'.join(','.join(str(s) for s in g) for g in haps.values()))
            
            w.writerow(row)
        
    if remove_singletons:
        logging.info("Number of singletons removed: {0}".format(num_singletons))

if __name__ == "__main__":
    def add_args(parser):
        input_file = parser.add_mutually_exclusive_group(required=True)
        input_file.add_argument("--infile", type="readable_file", metavar="FILE", default=None,
            help="File contaning genotypes.")
        input_file.add_argument("--bfile", metavar="PATH", default=None,
            help="PLINK binary file prefix.")
        parser.add_argument("--outfile", type="writeable_file", metavar="FILE",
            help="Output file.")
        parser.add_argument("--columns", metavar="INDEX", type='int_list', default=[0,1,2,3],
            help="List of four column indicies (0-based): snpID, chromosome, position and first " \
                 "genotype. Only used with --infile. (default=(0,1,2,3))")
        parser.add_argument("--delim", metavar="CHAR", default=",",
            help="Field delimiter in input and output files. (default=',')")
        parser.add_argument("--chrm", metavar="CHRM", default=None,
            help="Chromosome to process. Required if there is more than one chromosome in the input file.")
        parser.add_argument("--families", type="str_list", metavar="LIST", default=None,
            help="Comma-delimited list of family IDs to keep (only used with --bfile).")
        parser.add_argument("--haplotype-method", metavar="METHOD", choices=["c","a","l"], default="a",
            help="Method for identifying unique haplotypes. c=only collapse exact duplicates (most "\
                 "conservative), a=convert Ns to non-N calls (arbitrary), l=convert non-N calls to "\
                 "Ns (most liberal). (default=a)")
        parser.add_argument("--het", type="str_list", metavar="LIST", default=("H",),
            help="List of characters to treat as heterozygous.")
        parser.add_argument("--ignore-file", type="readable_file", metavar="FILE", default=None,
            help="File containing regions to ignore.")
        parser.add_argument("--interval-type", metavar="TYPE", 
            choices=["left", "right", "cores", "uber", "maxk"], default="maxk", 
            help="Algorithm for determining intervals. (default=maxk)")
        parser.add_argument("--missing-chars", type="str_list", metavar="LIST", default=("N","V","D"),
            help="List of characters to treat as missing information. The first character should "\
                 "be the 'no-call' character. If using a plink file (-b), this will be set to (3,1).")
        parser.add_argument("--maf", type=float, metavar="MAF", default=None,
            help="Minimum minor allele frequency of SNPs to include in analysis.")
        parser.add_argument("--max-n", metavar="FRACTION", type=float, default=0.1,
            help="Maximum fraction of N calls for a SNP. (default=0.2)")
        parser.add_argument("--max-snps", type="int_fmt", metavar="SNPS", default=None,
            help="Maximum number of SNPs that an interval can span.")
        parser.add_argument("--max-size", type="int_fmt", metavar="BP", default=None,
            help="Maximum size of an interval in bp.")
        parser.add_argument("--max-spacing", type="int_fmt", metavar="BP", default=None,
            help="Maximum distance between consecutive SNPs in the same interval.")
        parser.add_argument("--collapse-missing", action="store_true", default=False,
            help="Collapse haplotypes that differ only by missing characters. The default is to "\
                 "collapse only haplotypes that differ by N-calls.")
        parser.add_argument("--remove-singletons", action="store_true", default=False,
            help="Remove single-SNP intervals from output.")
        parser.add_argument("--write-haplotypes", action="store_true", default=False,
            help="Write unique haplotypes for each interval to the output file.")
        parser.add_argument("--write-trees", action="store_true", default=False,
            help="Compute and output a distance tree for each interval.")
            
    ns = parse(add_args)
    
    logging.debug("Parsing input file")
    
    if ns.infile is not None:
        reader = InfileReader(ns.infile, ns.chrm, ns.delim, ns.columns)
        het = ns.het
        missing = ns.missing_chars
    else:
        reader = PlinkReader(ns.bfile, ns.chrm, ns.families)
        het = ('1',)
        missing = ('3',)
    maf = ns.maf
    if maf is not None and maf <= 0:
        maf = None
    chrm = parse_input(reader, maf, ns.max_n, ns.ignore_file, het, missing)
    
    max_snps = ns.max_snps
    num_snps = len(chrm.snps)
    if max_snps is None:
        max_snps = num_snps
    else:
        assert max_snps <= num_snps
    
    max_size = ns.max_size
    chrm_size = chrm.snps[-1].position - chrm.snps[0].position + 1
    if max_size is None:
        max_size = chrm_size
    else:
        assert ns.max_size <= chrm_size
    
    max_spacing = ns.max_spacing
    if max_spacing is None:
        max_spacing = max_size - 1
    else:
        assert max_spacing < max_size
    
    logging.debug("Computing intervals")
    
    invs = chrm.get_intervals(ns.interval_type, 
        max_snps=max_snps, max_size=max_size, max_spacing=max_spacing)
    
    logging.info("Number of intervals: {0}".format(len(invs)))
    
    logging.debug("Writing output")
    
    if ns.write_trees:
        # This is a hack to avoid warnings from PyCogent about not having MPI enabled
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from cogent.phylo import nj
    
    if ns.collapse_missing:
        collapse_chars = missing + het
    else:
        collapse_chars = (missing[0],)
    
    write_output(chrm, invs, ns.outfile, ns.remove_singletons, ns.haplotype_method, 
        collapse_chars, ns.write_haplotypes, ns.write_trees)

class MockReader(object):
    def __init__(self, rows, header=None):
        self.rows = rows
        if header is None:
            header = (str(x) for x in xrange(0, len(rows[0][1])))
        self.samples = index_map(header)
        
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