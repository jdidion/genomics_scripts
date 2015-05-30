#! /usr/bin/env python

# Use a sliding window to classify segments of the genome as homo/heterozygous.
# Each window is identified by its call (0=heterozygous, 1=homozygous), starting
# position and percent homozygosity. A Smoother is then used to smooth the call
# sequence by changing small regions of fluctuation to the most appropriate call.
#
# The input CSV file has three header rows: SNP ID, chromosome position. Each
# subsequent row contains the sample ID followed by the genotypes for that 
# sample. There is also an optional samples file for mapping the sample ID in
# the genotype file to a sample name. This is a two column CSV with header
# ID,Name.
#
# Break each chromosome into intervals for a set of samples. An interval is the 
# smallest section of the chromosome that is either all homozygous or all 
# heterozygous for at least one sample. For example, consider the following two 
# samples:
#
# A: 00000011111111100000
# B: 11100000000111100000
#
# There are five intervals: 1-3, 4-6, 7-11, 12-15, and 16-20. Obviously, the
# number of intervals grows with the number of samples (or at least the number
# of haploytypes in the population).

from bitarray import bitarray
from csv import DictReader, reader, writer
import logging as log
import os
import re
import sys

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
from util.cl import parse
from util.collections import View, Range, zipiter, wrapiter
from util.gnet import karyotype
from util.io import linecount

class SNP(object):
    def __init__(self, name, chromosome, position, A, B):
        self.name = name
        self.chromosome = chromosome
        self.position = int(position)

class Chromosome(Range):
    def __init__(self, name, size, start, end):
        Range.__init__(self, start, end)
        self.name = name
        self.size = size

class Genotypes(object):
    def __init__(self, geno, fmt):
        ng = len(geno)
        self.H = bitarray(ng)
        self.N = bitarray(ng)
    
        for i in xrange(0, ng):
            # convert genotypes to two bitarrays: A) 0=N, 1=not N, B) 0=Het, 1=Hom
            g = int(geno[i])
            if fmt == 'num' and g > 0:
                g = 0 if g == 2 else 1
            self.H[i] = g
            self.N[i] = g+1
        
    def homozygosity(self, start, end):
        n = self.N[start:end]
        hom = float((n & self.H[start:end]).count(True))
        tot = n.count(True)
        return hom / tot

class Window(Range):
    def __init__(self, rng):
        Range.__init__(self, rng.start, rng.end)
        self.percs = []
    
    def append(self, perc):
        self.percs.append(perc)
    
    def __repr__(self):
        return "Window<{0}-{1}>".format(self.start, self.end)
    
class Region(object):
    def __init__(self, nsamples, start_pos):
        self.start_pos = start_pos
        self.end_pos = None
        self.calls = [] # 0 = het, 1 = hom
        for i in xrange(0, nsamples):
            self.calls.append(bitarray())
        self.windows = []
    
    def add_window(self, pos):
        w = Window(pos)
        self.windows.append(w)
        return w
    
    def append(self, sample, call):
        self.calls[sample].append(call)
    
    def start(self):
        if self.windows:
            return self.windows[0].start
        else:
            return None
    
    def end(self):
        if self.windows:
            return self.windows[-1].end
        else:
            return None
    
    def inner_range(self, incl=False):
        if self.windows:
            return Range(self.start(), self.end(), incl)
        else:
            return None
        
    def outer_range(self, incl=False):
        if self.start_pos and self.end_pos:
            return Range(self.start_pos, self.end_pos, incl)
        else:
            return None
            
    def __repr__(self):
        if self.windows:
            return "Region<{0}-{1}>".format(self.start(), self.end())
        else:
            return "Region<Empty>"

class Smoother(object):
    """
    Smooth out a sequence of 0's and 1's by requiring a minimum size run of either.
    TODO: it would be cool if I could make this operate directly on the bitarrays 
    """
    def __init__(self, window=20, max_iter=100):
        self.window = window
        self.max_iter = max_iter
        
        self.beg_exprs = (
            re.compile("^(0{1,%d})(1+)(.*)" % window),
            re.compile("^(1{1,%d})(0+)(.*)" % window))
        self.end_exprs = (
            re.compile("(.*)(1+)(0{1,%d})$" % window),
            re.compile("(.*)(0+)(1{1,%d})$" % window))
        self.mid_exprs = (
            re.compile("(.*1+)(0{1,%d})(1+.*)" % window),
            re.compile("(.*0+)(1{1,%d})(0+.*)" % window))

    def smooth(self, seq):
        if len(seq) < self.window:
            # Cannot smooth a sequence smaller than the smooth window
            if log.is_trace(): log.trace("Not smoothing sequence; too small")
            return seq
        
        def smooth_beg(seq):
            for e in self.beg_exprs:
                m = e.match(seq)
                if m:
                    a, b, rest = m.group(1, 2, 3)
                    if len(a) > len(b):
                        seq = "{0}{1}{2}".format(a, a[0] * len(b), rest)
                    else:
                        seq = "{0}{1}{2}".format(b[0] * len(a), b, rest)
                    return seq
            return None
        
        def smooth_end(seq):
            for e in self.end_exprs:
                m = e.match(seq)
                if m:
                    rest, a, b = m.group(1, 2, 3)
                    if len(a) > len(b):
                        seq = "{0}{1}{2}".format(rest, a, a[0] * len(b))
                    else:
                        seq = "{0}{1}{2}".format(rest, b[0] * len(a), b)
                    return seq
            return None
        
        def smooth_mid(seq):
            for e in self.mid_exprs:
                m = e.match(seq)
                if m:
                    beg, mid, end = m.group(1, 2, 3)
                    return "{0}{1}{2}".format(beg, end[0] * len(mid), end)
            return None
        
        for func in (smooth_beg, smooth_end, smooth_mid):
            for i in xrange(0, self.max_iter):
                newseq = func(seq)
                if newseq:
                    seq = newseq
                else:
                    break
        
        return seq
    
def read_samples(fname, delim="\t"):
    samples = {}
    if fname is not None:
        with open(fname, 'rU') as f:
            for d in DictReader(f, delimiter=delim):
                samples[d['ID']] = d['Name']
    return samples

def read_exclude(fname, delim="\t"):
    exclude = {}
    if fname is not None:
        with open(fname, 'rU') as f:
            for d in DictReader(f, delimiter=delim):
                c = exclude.setdefault(d['Chromosome'], [])
                c.append(Range(int(d['Start']), int(d['End']), True))
    return exclude

def classify_windows(snps, genos, exclude, window_size, slide, min_hom, smoothing_size):
    nsnp = len(snps)
    if log.is_debug(): log.debug("{0} snps".format(nsnp))

    nsam = len(genos)
    start = 0
    end = start + window_size - 1
    excl = exclude.next()
    
    regions = []
    region = None
    start_pos = snps[0].position

    while end < nsnp:
        pos = Range(snps[start].position, snps[end].position, True)
        if log.is_trace(): log.trace("Current region: {pos.start}-{pos.end}".format(pos=pos))
        
        if excl and excl.intersects(pos):
            if log.is_trace(): log.trace("Current region overlaps exclude region: {0}".format(excl))
            
            if region:
                region.end_pos = excl.start - 1
                if log.is_trace(): log.trace(
                    "Ending region {0} with {1} windows".format(region, len(region.windows)))
                region = None

            while end <= nsnp and excl >= snps[start].position:
                start += 1
            
            start_pos = excl.end + 1
            end = start + window_size - 1
            excl = exclude.next()
        
        else:
            if not region:
                region = Region(nsam, start_pos)
                regions.append(region)
            
            window = region.add_window(pos)
            if log.is_trace(): log.trace("Adding window for range {0}".format(pos))
            
            for i,g in enumerate(genos):
                hom = g.homozygosity(start, end)
                call = hom >= min_hom
                if log.is_trace(): log.trace(
                    "{0} markers are homozygous; call = {1}".format(hom, call))
                
                window.append(hom)
                region.append(i, call)
            
            start += slide
            end += slide

    if region:
        region.end_pos = snps[-1].position
        if log.is_trace(): log.trace(
            "Ending region {0} with {1} windows".format(region, len(region.windows)))
    
    if log.is_debug(): log.debug("Regions: {0}".format(regions))

    return regions

int_expr = re.compile("(0+|1+)")
def partition(r, ordered_samples, smoother, int_writer, min_writer):
    starts = set((0,))
    
    for i,c in enumerate(r.calls):
        sam = ordered_samples[i]
        seq = smoother.smooth(c.to01())
        ivs = int_expr.findall(seq)

        start = 0
        for j,iv in enumerate(ivs):
            call = iv[0]
            ilen = len(iv)
            next_start = start + ilen
            starts.add(next_start)
            
            # Each interval extends from the start position of the first window to the
            # to immediately before the start of the next window, unless it is the
            # first or last interval, in which case it starts at the region start/ends
            # at the region end
            int_writer.writerow((sam, 
                r.start_pos if j == 0 else r.windows[start].start,
                r.end_pos if j == len(ivs)-1 else r.windows[next_start].start - 1,
                ilen, call))
            
            start = next_start
    
    # consecutive interval starts across all samples define our minimal intervals
    starts = sorted(starts)
    lens = tuple(starts[i+1]-starts[i] for i in xrange(0, len(starts)-1))
    for i,iv in enumerate(zip(starts, lens)):
        start_pos = r.start_pos if i == 0 else r.windows[iv[0]].start
        end_pos = r.end_pos if i == len(lens)-1 else r.windows[starts[i+1]].start - 1
        min_writer.writerow([start_pos, end_pos, iv[1]] + [int(c[iv[0]]) for c in r.calls])
        
def main(argv=None):
    def add_opts(parser):
        parser.add_argument('-H', '--homozygosity_cutoff', type=int, default=97,
            help="Percent homozygosity required to declare a region homozygous.")
        parser.add_argument('-k', '--karyotype_file', type='readable_file', default=None,
            help="File containing the karyotype.")
        parser.add_argument('-m', '--smoothing_size', type=int, default=20,
            help="Window size to use when smoothing regions.")
        parser.add_argument('-s', '--sample_file', type='readable_file', default=None,
            help="File containing sample information to ")
        parser.add_argument('-w', '--window_size', type=int, default=300,
            help="Number of markers in sliding classification window.")
        parser.add_argument('-W', '--window_slide', type=int, default=1,
            help="Number of markers to slide the window.")
        parser.add_argument('-x', '--exclude_file', type='readable_file', default=None,
            help="File containing regions to ignore.")
        parser.add_argument('--genotype_format', choices=['num', 'bin'], default='num',
            help="Format of genotypes: num = -1/1/2/3/4, bin = -1=N, 0=hom, 1=het.")
        parser.add_argument('--max_smoothing_iterations', type=int, default=100,
            help="Maximimum iterations to spend smoothing a sequence.")
        parser.add_argument('genotype_file', type='readable_file',
            help="File containing genotypes for samples, one sample per row.")
        parser.add_argument('output_dir', type='writeable_dir',
            help="Directory to write results, one file per chromosome.")

    ns = parse(add_opts, args=argv)
    if log.is_debug(): log.debug("find_intervals.py called with args: {0}".format(ns))
    
    samples = read_samples(ns.sample_file)
    if log.is_debug(): log.debug("{0} samples".format(len(samples)))
    
    exclude = read_exclude(ns.exclude_file)
    if log.is_debug():
        if exclude:
            log.debug("Exclude regions on chromosomes {0}".format(exclude.keys()))
        else:
            log.debug("No exclude regions")
        
    chrom_sizes = dict((c.name, c.size) for c in karyotype(ns.karyotype_file))
    chromosomes = []
    ordered_samples = []
    genotypes = {}
    
    with open(ns.genotype_file, 'rU') as f:
        r = reader(f)
        
        # create list of tuples: SNPID, chromosome, position
        head_iter = zipiter(r.next() for i in xrange(0, 3))
        head_iter.next() # remove header column
        snps = map(lambda x: SNP(*x), head_iter)
        nsnp = len(snps)
        if log.is_debug(): log.debug("{0} snps".format(nsnp))

        chrom = None
        start = None
        for i in xrange(0, nsnp):
            s = snps[i]
            if s.chromosome != chrom:
                if start:
                    chromosomes.append(Chromosome(chrom, chrom_sizes[chrom], start, i))
                chrom = s.chromosome
                start = i
        chromosomes.append(Chromosome(chrom, chrom_sizes[chrom], start, nsnp))
        if log.is_debug(): log.debug("{0} chromosomes".format(len(chromosomes)))
        
        for sample in r:
            name = sample.pop(0) # pop off header column
            ordered_samples.append(samples.get(name, name))
            
            for chrom in chromosomes:
                genotypes.setdefault(chrom.name, []).append(
                    Genotypes(chrom.slice(sample), ns.genotype_format))
                
    min_hom = float(ns.homozygosity_cutoff) / 100
    smoother = Smoother(ns.smoothing_size, ns.max_smoothing_iterations)
    
    if not os.path.exists(ns.output_dir):
        os.makedirs(ns.output_dir)
    
    for chrom in chromosomes:
        log.info("Processing chromosome {0}".format(chrom.name))

        # break genotypes into regions, scan each region using a sliding window
        # and classify each window as homozygous or heterozygous
        regions = classify_windows(View(snps, chrom), genotypes[chrom.name], 
            wrapiter(exclude.get(chrom.name, None)), 
            ns.window_size, ns.window_slide, min_hom, ns.smoothing_size)
        
        int_file = os.path.join(ns.output_dir, "intervals_chr%s.csv" % chrom.name)
        min_file = os.path.join(ns.output_dir, "minimal_chr%s.csv" % chrom.name)
        
        # smooth out each region and partition into hom/het intervals
        with open(int_file, 'w') as iout, open(min_file, 'w') as mout:
            int_writer = writer(iout)
            int_writer.writerow(('Sample','Start','End','Martkers','Call'))
            
            min_writer = writer(mout)
            min_writer.writerow(['Start','End','Markers'] + ordered_samples)
            
            for r in regions:
                partition(r, ordered_samples, smoother, int_writer, min_writer)
                
if __name__ == '__main__':
    main()