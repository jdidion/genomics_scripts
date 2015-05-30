# Utility commands for working with Plink
# This library depends on two module variables:
# * PLINK_CMD (path to plink version >= 1.9)
# * PLINK_OLD (path to plink version 1.07)
from bisect import bisect_left, bisect_right
from collections import OrderedDict, defaultdict, namedtuple
from csv import reader, writer
from itertools import izip
import logging
from math import floor,ceil
import os
from random import choice

import numpy as np
import scipy.stats as st

from util.collections import which
from util.io import FileCloser, csv_to_dict, csv_to_table, safe_read_file, safe_read_file_array, mkdir, filename
from util.misc import cumsum

PLINK_CMD = None
PLINK_OLD = None

# TODO: this function should move to the conversion library
# TODO: need to rename samples since plink files are space-delimited
def convert_raw_to_transposed(geno_file, sample_file, snps_file, plink_tped_file, plink_tfam_file, 
        chromosomes, vinos_as_synthetic):
    """
    Convert raw input to Plink transposed format. This method expects three files:
    1. genotypes.csv: The genotypes file produced by MouseDivGeno. The first row contains the sample 
       names. For each successive row, the first column is the SNP ID, and the following columns are 
       genotypes in numeric format (1=A, 2=H, 3=B, 4=V, 5=N).
    2. snps.csv: Each row has six columns: SNP ID (must match the SNP ID in the genotype file), 
       chromosome, physical position, genetic position (in cM), reference allele, variant allele.
    3. samples.csv: Each row has at least 7 columns: sample ID, sex, population/family, mother ID, 
       father ID, 2N, translocations. Sex should be coded as 0=unknown, 1=male, 2=female. Population 
       is useful even if it is just a guess. An example for wild mice would be to consider each 
       subspecies a population. If you don't want to provide this, just give every sample the same
       population id. Parent IDs should be 0 if you don't want to specify parents. 2N is the
       diploid number and should be 40 unless known to be different. If 2N < 40, the translocation 
       column should specific translocations as A.B, where A and B are the fused chromosomes. The list 
       of translocations should be space-delimited. Heterozygous translocations are denoted by an 
       asterisk (e.g. A.B*). Any columnts beyond 7 are treated as phenotype values.
    """
    
    samples = csv_to_table(sample_file)
    sample_header = samples[0]
    assert len(sample_header) >= 5, "Too few columns in sample file"
    samples = samples[1:]
    nsamples = len(samples)
    sample_rng = range(nsamples)

    chrms = set(map(str, chromosomes))
    bp = 1/1000000 # genetic distance covered by 1 bp

    with open(geno_file, "rU") as gfile, open(snps_file, "rU") as sfile, open(plink_tped_file, "w") as tfile:
        geno_reader = reader(gfile)
        geno_header = geno_reader.next()[1:]
        assert len(geno_header) >= nsamples, "Too few columns in genotype file"
        ixs = tuple(which(s[0], geno_header) for s in samples)

        snp_reader = reader(sfile)
        snps_header = snp_reader.next()
        assert len(snps_header) == 6, "Incorrect number of columns in snps file"

        tped_writer = writer(tfile, delimiter=" ")

        genos = [0] * (nsamples * 2)
        vinos = [0] * (nsamples * 2)
        # TODO: we should probably instead assign an allele pair at random using
        # probabilities from the standard transition/transversion distribution
        V = ("C", "C", "C", "G", "0")
        
        for row,data in enumerate(izip(geno_reader, snp_reader)):
            sid, chrm, ppos, gpos, ref, var = data[1]
            if chrm not in chrms: continue
            assert sid == data[0][0], "SNP does not match genotype at row {0}".format(row)
            geno = data[0][1:]
            
            A = (ref,ref,var,"0","0")
            B = (ref,var,var,"0","0")
            vino = False
            for i,x in enumerate(ixs):
                g = int(geno[x]) - 1
                genos[(i*2)] = A[g]
                genos[(i*2)+1] = B[g]
                if vinos_as_synthetic:
                    vinos[(i*2)] = V[g]
                    vinos[(i*2)+1] = V[g]
                    if g == 3: 
                        vino = True

            tped_writer.writerow([chrm, sid, gpos, ppos] + genos)
            if vino:
                tped_writer.writerow([chrm, "{0}V".format(sid), float(gpos) + bp, int(ppos) + 1] + vinos)

    with open(plink_tfam_file, "w") as tfile:
        tfam_writer = writer(tfile, delimiter=" ")
        for i in sample_rng:
            sample = samples[i]
            tfam_writer.writerow((sample[2], sample[0], sample[3], sample[4], sample[1], "-9"))

## Command formatters for single-thread operations

def convert_transposed_to_binary(infile, outfile=None):
    if outfile is None:
        outfile = filename(infile)
    return "{cmd} --tfile {infile} --out {outfile} --make-bed".format(cmd=PLINK_CMD, infile=infile, outfile=outfile)

def filter_genotypes(infile, outfile, filters):
    return ("{0} --bfile {1} --out {2} --make-bed".format(PLINK_CMD, infile, outfile),) + tuple(
        "--{0} {1}".format(k, v) for k,v in args.filters.iteritems())

def format_pairwise_ld_command(bfile, outfile, ld_window_r2=0, sample_file=None, snp_file=None, 
        inter_chr=True, apply_filters=False, min_maf=0.01, max_uncalled=0.2):
    """Create a plink command to generate an LD matrix for all possible SNP comparisons."""

    cmd = "{0} --bfile {1} --out {2} --r2 --matrix".format(PLINK_CMD, bfile, outfile)
    if sample_file:
        cmd += " --keep {0}".format(sample_file)
    if snp_file:
        cmd += " --ld-snp-list {0}".format(snp_file)
    if apply_filters:
        cmd += " --maf {0} --geno {1}".format(min_maf, max_uncalled)
    return cmd
            
## Command iterators for threadable operations

def generate_stats(infile, outfile, stats):
    if outfile is None:
        outfile = filename(infile)
    av = dict(cmd=PLINK_CMD, infile=infile, outfile=outfile)
    cmd = "{cmd} --bfile {infile} --out {outfile} --{stat}"
    for s in args.stats:
        av["stat"] = s
        yield(av, cmd.format(**av), {})

FORMAT = dict(bed="make-bed", ped="recode", tped="transpose")

def split(infile, outdir, chromosomes=range(1,20)+["X","Y","M"], output_format="bed"):
    av = dict(infile=infile)
    outfile = filename(infile)
    for chrm in chromosomes:
        chrm_dir = os.path.join(outdir, "chr{0}".format(chrm))
        mkdir(chrm_dir, overwrite=False)
        av["cmd"] = PLINK_CMD
        av["outfile"] = os.path.join(chrm_dir, outfile)
        av["chrm"] = chrm
        av["format"] = FORMAT[output_format]
        cmd = "{cmd} --bfile {infile} --chr {chrm} --out {outfile} --{format}".format(**av)
        yield (av, cmd, {})

def per_chrm_iter(infile, outdir, fn, chromosomes=range(1,20)+["X","Y","M"], *args, **kwargs):
    fname = filename(infile)
    for chrm in chromosomes:
        chrm_dir = os.path.join(outdir, "chr{0}".format(chrm))
        chrm_file = os.path.join(chrm_dir, fname)
        cmd = fn(chrm_file, chrm_file, *args, **kwargs)
        yield (dict(infile=chrm_file, outfile=chrm_file, chrm=chrm), cmd, {})

## Calculate genome-wide pairwise LD matrix.
# These methods support a framework for doing LD computation in parallel and merging the results.
# First, SNPs are partitioned into windows of a certain size. Then plink is used to calculate
# pairwise LD between all SNPs in a window and every other SNP. The resulting files (in matrix
# format) are processed and summary files are generated. Plots can be created from the summaries.
# Files are expected to have specific names of the form: <outdir>/chr<chrm>/win_<window>.<ext>.

def get_window(x, size): 
    return floor(x / size) * size

# TODO: take a chromosomes paramater and only make windows for those chromosomes
def partition_snps(bimfile, outdir, window_size=500000):
    """
    Read a plink .bim file and partition SNPs into groups defined by window_size. This is to
    faciliate parallization of plink LD calculation. Also writes two summary csv files:
    1) chromosomes: chrm, start_win_idx, end_win_idx, start_snp_idx, end_snp_idx
    2) windows: chrm, window, start_snp_idx, end_snp_idx.
    """
    
    chr_file = os.path.join(outdir, "chromosomes.csv")
    win_file = os.path.join(outdir, "windows.csv")
    pos_file = os.path.join(outdir, "positions.txt")
    
    with open(bimfile, "rU") as b, open(chr_file, "w") as c, open(win_file, "w") as w, open(pos_file, "w") as p:
        chr_writer = writer(c)
        win_writer = writer(w)

        cur_chrm = None
        cur_chrdir = None
        cur_win = None
        cur_out = None
        
        prev_pos = 0
        cur_win_idx = 0
        chr_start_win_idx = 0
        chr_start_snp_idx = 0
        win_start_snp_idx = 0
        
        for i,line in enumerate(reader(b, delimiter="\t")):
            chrm = line[0]
            new_chrm = False
            if i > 0 and chrm != cur_chrm:
                chr_writer.writerow((cur_chrm, chr_start_win_idx, cur_win_idx, chr_start_snp_idx, i-1))
                new_chrm = True
                chr_start_snp_idx = i
                chr_start_win_idx = cur_win_idx + 1
                    
            pos = int(line[3])
            if new_chrm:
                prev_pos = 0
            elif pos < prev_pos:
                raise Exception("bimfile is not sorted")
            else:
                prev_pos = pos
            p.write(str(pos))
            p.write("\n")
            
            win = int(get_window(pos, window_size))
            new_win = False
            if new_chrm or win != cur_win:
                if cur_win is not None:
                    win_writer.writerow((cur_chrm, cur_win, win_start_snp_idx, i-1))
                new_win = True
                cur_win_idx += 1
                win_start_snp_idx = i
            
            if new_chrm or new_win:
                if cur_out is not None:
                    cur_out.close()

                if new_chrm or i == 0:
                    cur_chrdir = os.path.join(outdir, "chr{0}".format(chrm))
                    mkdir(cur_chrdir, overwrite=True)
                    cur_chrm = chrm

                outfile = os.path.join(cur_chrdir, "win_{0}.snps".format(win))
                cur_out = open(outfile, "w")
                cur_win = win
            
            snp = line[1]
            cur_out.write(snp)
            cur_out.write("\n")
        
        chr_writer.writerow((cur_chrm, chr_start_win_idx, cur_win_idx, chr_start_snp_idx, i))
        win_writer.writerow((cur_chrm, cur_win, win_start_snp_idx, i))
        cur_out.close()
        
    return (chr_file, win_file)

def pairwise_ld_command_iter(bfile, outdir, window_file, **kwargs):
    """Create a command iterator over a set of SNP files (created by partition_snps)."""
    
    with open(window_file, "rU") as w:
        for row in reader(w):
            chr_dir = os.path.join(outdir, "chr{0}".format(row[0]))
            snp_file = os.path.join(chr_dir, "win_{0}.snps".format(row[1]))
            out_file = os.path.join(chr_dir, "win_{0}".format(row[1]))
            yield(dict(chrm=row[0], window=row[1]), 
                format_pairwise_ld_command(bfile, out_file, snp_file=snp_file, **kwargs), {})
        
class Bin(object):
    def __init__(self, r2_bin_size):
        self.count = 0
        self.sum_r2 = 0
        self.max_r2 = 0
        # We calculate percentiles using an approximation based on binning values.
        # The average error of this method is +/- (r2_bin_size / 2)
        self.r2_bin_size = r2_bin_size
        self.bin_r2 = [0] * int((1 / r2_bin_size) + 1)

    def add(self, r2):
        self.count += 1
        self.sum_r2 += r2
        if self.max_r2 < r2:
            self.max_r2 = r2
        self.bin_r2[int(floor(r2 / self.r2_bin_size))] += 1
        
    def as_list(self, percentile=95, prec=6):
        bin_avg = round(self.sum_r2 / self.count, prec)
        bin_max = round(self.max_r2, prec)
        bin_pct = round(self._percentile(percentile), prec)
        return [self.count, bin_avg, bin_max, bin_pct]
    
    def _percentile(self, pct):
        """Estimate the value at the given percentile from binned data.
        """
        c = tuple(cumsum(self.bin_r2))
        x = c[-1] * pct / 100.0
        l = bisect_left(c, x)
        if l == 0:
            b = x / c[l]
        elif l == len(c) - 1:
            b = l
        else:
            b = l + ((x - c[l - 1]) / (c[l] - c[l - 1]))
        return b * self.r2_bin_size

def process_ld_files(outdir, percentile=95, interval_bin_size=1000, r2_bin_size=0.01, unlinked=True):
    """
    Process plink.ld matrix files and generate summary files. The files created are:
    1) NxN matrix of mean, max and Xth percentile r^2 values for each pairwise window comparison, 
       where N is the number of windows and X is specified by the percentile argument. To speed 
       computation, only the upper right triangle of each matrix is computed.
    2) The mean, max and Xth percentile r^2 values for each window. The summaries are computed for
       both the whole genome and for only those windows that are not on the same chromosome as the
       window.
    3) The count, mean, max and Xth percentile r^2 value for bins of inter-SNP distance (determined 
       by the interval_bin_size argument).
    4) The count of r^2 values in bins specified by r2_bin_size.
    """
    
    chr_file = os.path.join(outdir, "chromosomes.csv")
    chromosomes = csv_to_dict(chr_file, convert=lambda row: map(int, row))
    nchr = len(chromosomes)
    
    win_file = os.path.join(outdir, "windows.csv")
    windows = csv_to_table(win_file, convert=lambda row: map(int, row))
    nwin = len(windows)
    
    pos_file = os.path.join(outdir, "positions.txt")
    positions = safe_read_file_array(pos_file, convert=int)
    nsnps = len(positions)
    
    stats = OrderedDict()
    stats["mean"] = np.mean
    stats["max"] = np.max
    stats["pct"] = lambda m: st.scoreatpercentile(m, percentile)
    
    def summarize(m):
        a = m.ravel()
        a = a[~np.isnan(a)]
        return [fn(a) for fn in stats.values()] if len(a) > 0 else ["nan"] * len(stats)

    def myopen(f, header, closer):
        fn = os.path.join(outdir, f)
        fp = open(fn, "w")
        wr = writer(fp)
        if header:
            wr.writerow(header)
        closer.add(fp)
        return wr
    
    mat_header = ['',] + map(lambda w: "{0}_{1}".format(w[0], w[1]), windows)
    win_header = ("chrm", "window", "snps", "all_mean", "all_max", "all_pct")
    if unlinked:
        win_header += ("unlinked_mean", "unlinked_max", "unlinked_pct")
    
    closer = FileCloser()
    mat_writers = {}
    for x in stats.keys():
        mat_writers[x] = myopen("{0}_r2_matrix.csv".format(x), mat_header, closer)
    win_writer = myopen("window_summary.csv", win_header, closer)
    
    inter_snp = defaultdict(lambda: Bin(r2_bin_size))
    
    r2prec = len(str(r2_bin_size)) - 2 # precision of r2 bins
    r2_count_all = defaultdict(int)
    r2_count_intra = defaultdict(int)
    r2_count_inter = defaultdict(int)
    
    # process ld files one window at a time
    for i,win1 in enumerate(windows):
        win1_name = "{0}_{1}".format(win1[0], win1[1])
        chrm, first_win, last_win, first_snp, last_snp = chromosomes[win1[0]]
        
        ld_file = os.path.join(outdir, "chr{0}".format(win1[0]), "win_{0}.ld".format(win1[1]))
        logging.debug("Processing {0}".format(ld_file))
        matrix = np.genfromtxt(ld_file)
        if len(matrix.shape) == 1:
            matrix = matrix.reshape(1, matrix.shape[0])
        logging.debug("Matrix shape: {0}".format(matrix.shape))
        
        # calculate whole-genome stats
        win_row = [win1[0], win1[1], win1[3] - win1[2] + 1]
        win_row.extend(summarize(matrix))
        
        # unlinked stats exclude the current chromosome
        if unlinked:
            before = matrix[..., 0:windows[first_win][2]] if first_win > 0 else None
            after = matrix[..., windows[last_win + 1][2]:nsnps] if last_win + 1 < nwin else None
            if before is None:
                unlinked = after
            elif after is None:
                unlinked = before
            else:
                unlinked = np.concatenate((before, after), axis=1)
            win_row.extend(summarize(unlinked))
        
        win_writer.writerow(win_row)
        
        # calculate pairwise window stats
        mat_rows = {}
        for x in stats.keys():
            mat_rows[x] = [win1_name] + ([0] * i)
            
        for win2 in windows[i:nwin]:
            s = summarize(matrix[..., win2[2]:(win2[3] + 1)])
            for i in xrange(0, len(stats)):
                v = s[i] if s else "nan"
                mat_rows[stats.keys()[i]].append(v)
        
        for x in stats.keys():
            mat_writers[x].writerow(mat_rows[x])
        
        # add to the inter-SNP distance summary
        for row in xrange(0, matrix.shape[0]):
            snp_idx = win1[2] + row
            pos = positions[snp_idx]
            for chr_idx in xrange(snp_idx + 1, last_snp + 1):
                r2 = matrix[row, chr_idx]
                if not np.isnan(r2):
                    dist = int(get_window(positions[chr_idx] - pos, interval_bin_size))
                    inter_snp[dist].add(r2)
        
        # add to the r^2 histogram
        for c2 in chromosomes.values():
            chrm2 = c2[0]
            for r2 in matrix[..., c2[1]:(c2[2] + 1)].ravel():
                if not np.isnan(r2):
                    r2_win = round(get_window(r2, r2_bin_size), r2prec)
                    r2_count_all[r2_win] += 1
                    if chrm == chrm2:
                        r2_count_intra[r2_win] += 1
                    else:
                        r2_count_inter[r2_win] += 1
        
        del matrix
    
    # Write inter-SNP and r2 count files
    int_header = ("size", "count", "mean", "max", "pct{0}".format(percentile))
    w = myopen("inter_snp.csv", int_header, closer)
    for b in sorted(inter_snp.keys()):
        w.writerow([b] + inter_snp[b].as_list(percentile))
    
    r2_header  = ("r2", "count")
    def write_r2_hist(d, f):
        w = myopen(f, r2_header, closer)
        for r in sorted(d.keys()):
            w.writerow((r, d[r]))
    
    write_r2_hist(r2_count_all, "r2_hist_all.csv")
    write_r2_hist(r2_count_inter, "r2_hist_inter.csv")
    write_r2_hist(r2_count_intra, "r2_hist_intra.csv")
    
    closer.close()

def process_local_ld_file(infile, outfile, percentile=95, bin_size=1000):
    # Process an LD file generated using the --r2 command with paramters to look for LD in cis.
    # SNP pairs are always assumed to be on the same chromosome. The r^2 for each SNP pair is
    # added to a bin, and then the Nth percentile, mean and SD for each bin are printed to outfile.
    bins = defaultdict(lambda: [])
    with open(infile, "rU") as i:
        r = reader(i, delimiter=" ", skipinitialspace=True)
        r.next() # skip header
        for row in r:
            p1 = int(row[1])
            p2 = int(row[4])
            d = abs(p2-p1)
            b = int(ceil(float(d) / bin_size) * bin_size)
            r2 = float(row[6])
            bins[b].append(r2)
    
    with open(outfile, "w") as o:
        w = writer(o)
        for b in sorted(bins.keys()):
            p = st.scoreatpercentile(bins[b], percentile)
            m = np.mean(bins[b])
            s = np.std(bins[b])
            w.writerow((b, p, m, s))

def format_allele_freq_command(infile, outfile=None):
    return "{0} --bfile {1} --freq --out {2}".format(PLINK_CMD, infile, outfile or infile)

def format_hap_freq_command(infile, outfile=None):
    return "{0} --bfile {1} --hap-window 2 --hap-freq --out {2}".format(PLINK_CMD, infile, outfile or infile)

def format_ld_blocks_command(infile, outfile=None, window_kb=500, apply_filters=False, min_maf=0.01, max_uncalled=0.2):
    cmd = "{0} --bfile {1} --out {2} --blocks --ld-window-kb {3}".format(PLINK_CMD, infile, outfile or infile, window_kb)
    if apply_filters:
        cmd += " --maf {0} --geno {1}".format(min_maf, max_uncalled)
    return cmd

def block_ld(line_iter, infile, outdir, digits=5):
    """
    Generate commands to compute r-squared values for all pairwise combinations 
    of SNPs in each block and the subsequent block.
    """

    cur_chrm = None
    cur_block = 1
    prev_snps = None
    all_snps = None
    summary = None
    chrm_dir = None
    
    def end_chrm():
        exclude_outfile = os.path.join(chrm_dir, "block_snps.txt")
        include_outfile = os.path.join(chrm_dir, "singleton_snps.txt")
        summary_file = os.path.join(chrm_dir, "block_ld_summary.txt")

        with open(summary_file, "w") as s:
            w = writer(s)
            w.writerow(("block1", "block2", "block1_snps", "block2_snps", "snp_file", "ld_file"))
            w.writerows(summary)
        with open(exclude_outfile, "w") as o:
            o.write("\n".join(all_snps))
        av = dict(infile=infile, chrm=cur_chrm, exclude=exclude_outfile, include=include_outfile)
        return [av, "{0} --bfile {1} --chr {2} --exclude {3} --write-snplist --out {4}".format(
            PLINK_CMD, infile, cur_chrm, exclude_outfile, include_outfile), {}]
    
    for row in reader(line_iter, delimiter=" ", skipinitialspace=True):
        # skip header
        if row[0] == "CHR":
            continue
            
        chrm, start, end, kb, nsnps, snps = row
        snps = snps.split("|")
        if cur_chrm != chrm:
            if cur_chrm is not None:
                yield end_chrm()
            
            chrm_dir = os.path.join(outdir, "chr{0}".format(chrm))
            mkdir(chrm_dir, overwrite=False)
            
            cur_chrm = chrm
            cur_block = 1
            all_snps = set(snps)
            summary = []
        else:
            fname = os.path.join(chrm_dir, "block_{0}".format(str(cur_block).zfill(digits)))
            outfile = "{0}.snps".format(fname)
            with open(outfile, "w") as o:
                o.write("\n".join(prev_snps + snps))
            summary.append((cur_block, cur_block+1, "|".join(prev_snps), "|".join(snps), outfile, "{0}.ld".format(fname)))
            
            av = dict(infile=infile, outfile=outfile, outname=fname)
            yield [av, "{0} --bfile {1} --extract {2} --r2 --matrix --out {3}".format(PLINK_CMD, infile, outfile, fname), {}]
            
            cur_block += 1
            all_snps.update(snps)
        
        prev_snps = snps

    end_chrm()

def prune_blocks(outdir, chromosomes=range(1,20)+["X","Y","M"], digits=5, prec=3):
    np.set_printoptions(precision=prec)
    snps_file = os.path.join(outdir, "ld_pruned_snp_list.txt")
    with open(snps_file, "w") as s:
        for chrm in chromosomes:
            chrm_dir = os.path.join(outdir, "chr{0}".format(chrm))
            if not os.path.isdir(chrm_dir):
                continue
                
            old_summary_file = os.path.join(chrm_dir, "block_ld_summary.txt")
            new_summary_file = os.path.join(chrm_dir, "block_ld_selected_snps.txt")
            
            with open(old_summary_file, "rU") as i, open(new_summary_file, "w") as o:
                summary = writer(o)
                summary.writerow(("block1", "block2", "block1_snps", "block2_snps", "block1_selected", "block1_r2", "block2_selected", "block2_r2"))
                min_snp = None
                
                r = reader(i)
                r.next() # skip header
                for row in r:
                    block1, block2, block1_snps, block2_snps, snp_file, ld_file = row
                    row = row[0:4]
                    
                    snps1 = block1_snps.split("|")
                    snps2 = block2_snps.split("|")
                    l = len(snps1)
                    matrix = np.genfromtxt(ld_file)[:l,l:]
                    
                    # find the SNP in block 1 with the lowest LD with block 2 SNPs
                    min_idx = np.argmin(np.apply_along_axis(np.sum, 1, matrix))
                    min_snp = snps1[min_idx]
                    row += [min_snp, np.round(np.mean(matrix[min_idx,:]), prec)]
                    
                    s.write(min_snp)
                    s.write("\n")
                    
                    # find the SNP in block 2 2ith the lowest LD with block 1 SNPs
                    min_idx = np.argmin(np.apply_along_axis(np.sum, 0, matrix))
                    min_snp = snps2[min_idx]
                    row += [min_snp, np.round(np.mean(matrix[:,min_idx]), prec)]
                    
                    summary.writerow(row)
                    
                # add snp from last block
                s.write(min_snp)
                s.write("\n")
                
                # append all snps that were not in blocks
                s.write(safe_read_file(os.path.join(chrm_dir, "singleton_snps.txt.snplist")))

# Before make_blocks, run commands such as the following:
# ~/sci/common/exe/util/exec.py -v chrm=F -f thread -k threads=14 --log_level DEBUG {PLINK_CMD} --bfile ../mm --chr "{chrm}" --r2 in-phase dprime with-freqs --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0 --freqx --out "chr{chrm}/mm"
# {plink_old_cmd} --bfile ../mm --chr "{chrm}" --hap-window 2 --hap-freq --out "chr{chrm}/mm"

def make_blocks(chrm, path, r2_lower=0.15, r2_upper=0.95, min_consecutive_unlinked=2, min_maf_4g=0.05):
    """Given .ld, .frq.hap and .frqx files for a chromosome, identify haplotype blocks and tag SNPs. 
    
    Starting with the first SNP (A), all pairwise comparisons between A and subsequent SNPs (B) are 
    examined. If all four diplotypes are observed with frequency greater than the specified 
    threshold (min_maf_4g) for any A-B, A and B must be separated by recombination due to the 
    four-gamete test and the block may maximally extend from A to the SNP prior to B. Next, if the 
    r-squared value is less than the lower threshold (r2_lower) for the specified number of 
    consecutive SNPs (min_consecutive_unlinked) then the provisional end of the block is the last 
    linked SNP. Any remaining SNPs in the block that have an r-squared value below the threshold
    are not considered in further analyses. Next, if there a pair is encountered with an r-squared 
    value that exceeds r2_upper, then B is considered the new tag SNP for the block. In that case, if 
    the block does not have a hard distal boundary (due to the four-gamete test), the block may be 
    extended by markers that are in LD to the new tag SNP.
    
    Two output files are written: a PLINK include file with the names of the tag SNPs, and 
    a more detailed block file with the start and end of each block, the tag SNP, and a flag column 
    that is a bitwise combination of the following: 
    
    1 = singleton block
    2 = block ends due to four-gamete test
    """

    class SNP(object):
        def __init__(self, name, missing):
            self.name = name
            self.missing = missing
            self.breakpoint = False
            self.pos = None
            self.maf = None
        
        @property
        def seen(self):
            return self.pos is not None
        
        def __str__(self):
            return "SNP<{0},{1}>".format(self.name, self.pos)
    
    # first read all SNPs into a dict
    logging.debug("Reading SNPs")
    snps = OrderedDict()
    snp_file = "{0}.frqx".format(path)
    with open(snp_file, "rU") as f:
        snp_reader = reader(f, delimiter="\t")
        header = snp_reader.next() # skip header
        for snp in snp_reader:
            name = snp[1]
            snps[name] = SNP(name, int(snp[9]))

    # second, look for any diplotypes that fail the four-gamete test
    logging.debug("Reading diplotypes")
    hap_file = "{0}.frq.hap".format(path)
    with open(hap_file, "rU") as h:
        hap_reader = reader(h, delimiter=" ", skipinitialspace=True)
        header = hap_reader.next() # skip header
        locus = None
        F = None
        snp_idx = -1
        for hap in hap_reader:
            if locus != hap[0]:
                locus = hap[0]
                F = []
                snp_idx += 1
            F.append(float(hap[2]))
            if len(F) == 4 and all(f > min_maf_4g for f in F):
                k = snps.keys()[snp_idx]
                snps[k].breakpoint = True
    
    
    ld_file = "{0}.ld".format(path)
    out_file = "{0}.blocks".format(path)
    inc_file = "{0}.tag_snps".format(path)
    
    with open(ld_file, "rU") as ld, open(out_file, "w") as out, open(inc_file, "w") as inc:
        ld_reader = reader(ld, delimiter=" ", skipinitialspace=True)
        ld_reader.next() # skip header
        
        out_writer = writer(out, delimiter="\t")
        out_writer.writerow(("Chr", "Start_Name", "Start_Pos", "End_Name", "End_Pos", "Tag_Name", "Tag_Pos", "Tag_Missing", "Flag"))
        
        class SNP_Pair(object):
            def __init__(self, row):
                self.A = snps[row[2]]
                if not self.A.seen:
                    self.A.pos = int(row[1])
                    self.A.maf = float(row[3])
                
                self.B = snps[row[6]]
                if not self.B.seen:
                    self.B.pos = int(row[5])
                    self.B.maf = float(row[8])
            
                self.r2 = float(row[9])
    
        class Block(object):
            def __init__(self, pair):
                self.snps = []
                self.cur_tag = None
                self.best_tag = None
                self.ended = self.fourg = False
                self.unlinked = 0
                
                self.add_snp(pair.A, True)
                if not self.ended:
                    self.next_snp(pair.B, pair.r2)
            
            def __str__(self):
                return "Block<Proximal={0},Distal={1},CurrentTag={2},BestTag={3},NumSnps={4},Size={5},Ended={6},4G={7}>".format(
                    self.start, self.end, self.cur_tag, self.best_tag, len(self.snps), self.end.pos-self.start.pos+1, self.ended, self.fourg)
            
            @property
            def singleton(self):
                return self.ended and self.size == 1
            
            @property
            def start(self):
                return self.snps[0]
            
            @property
            def end(self):
                return self.snps[-1]
            
            @property
            def size(self):
                return len(self.snps)
            
            def add_snp(self, snp, is_tag=False, perfect_linkage=False):
                self.snps.append(snp)
                if is_tag:
                    self.cur_tag = snp
                    if self.best_tag is None or (
                            self.cur_tag.missing < self.best_tag.missing) or (
                            self.cur_tag.missing == self.best_tag.missing and self.cur_tag.maf > self.best_tag.maf):
                        self.best_tag = snp
                    elif perfect_linkage:
                        self.best_tag = choice((snp,self.best_tag))
                    self.unlinked = 0
                if snp.breakpoint:
                    self.ended = self.fourg = True
            
            def next_snp(self, snp, r2):
                if r2 > r2_lower:
                    self.add_snp(snp, r2 > r2_upper, r2 == 1.0)
                    self.unlinked = 0
                else:
                    self.unlinked += 1
                    if self.unlinked == min_consecutive_unlinked:
                        self.ended = True
            
            def next_pair(self, pair):
                if not self.ended:
                    if pair.A == self.cur_tag:
                        self.next_snp(pair.B, pair.r2)
                    elif pair.A.pos > self.cur_tag.pos:
                        self.ended = True
                        
                if self.ended and pair.A.pos > self.end.pos:
                    self.write()
                    return Block(pair)
                else:
                    return self
            
            def write(self):
                start = self.start
                end = self.end
                tag = self.best_tag
                flags = (1 if self.singleton else 0) | (2 if self.fourg else 0)
                out_writer.writerow((chrm, start.name, start.pos, end.name, end.pos, tag.name, tag.pos, tag.missing, flags))
                inc.write(tag.name)
                inc.write("\n")
        
        logging.debug("Finding blocks")
        
        block = Block(SNP_Pair(ld_reader.next()))
        
        for row in ld_reader:
            block = block.next_pair(SNP_Pair(row))
        
        # write last block
        block.write()

def format_merge_tag_lists_command(outdir, name, chromosomes=range(1,20)+["X","Y","M"]):
    infiles = " ".join(os.path.join(outdir, "chr{0}".format(chrm), "{0}.tag_snps".format(name)) for chrm in chromosomes)
    outfile = os.path.join(outdir, "{0}.tag_snps".format(name))
    return "cat {0} > {1}".format(infiles, outfile)

# This produces an invalid format. Use PGDSpider instead.
def format_structure_command(bfile, outdir, name, chromosomes=range(1,20)+["X"]):
    av = dict(cmd=PLINK_CMD, bfile=bfile)
    for chrm in chromosomes:
        outfile = os.path.join(outdir, "{0}.chr{1}".format(name, chrm))
        av["chrm"] = chrm
        av["outfile"] = outfile
        cmd = "{cmd} --bfile {bfile} --chr {chrm} --recode structure --out {outfile}".format(**av)
        yield [av, cmd, {}]

# Currently, PLINK uses physical distance for the map row in structure files.
# This function replaces it with genetic distance.
def replace_structure_maps(bfile, outdir, name):
    map_file = "{0}.bim".format(bfile)
    
    def my_fmt(f):
        f = round(f,5)
        if f == 0:
            return "0.0"
        else:
            return format(f, 'f').rstrip('0')
    
    newlines = {}
    with open(map_file, "rU") as m:
        gmap = defaultdict(lambda: [])
        for row in reader(m, delimiter="\t"):
            gmap[row[0]].append(float(row[2]))
        for chrm, cm in gmap.iteritems():
            dif = ("-1",) + tuple(my_fmt(cm[i] - cm[i-1]) for i in range(1, len(cm)))
            newlines[chrm] = "\t\t{0}\n".format(" ".join(dif))
    
    for chrm, newline2 in newlines.iteritems():
        if chrm == "23":
            chrm = "X"
        structure_file = os.path.join(outdir, "{0}.chr{1}.recode.strct_in".format(name, chrm))

        with open(structure_file, "rU") as s:
            line1 = s.readline()
            oldline2 = s.readline()
            rest = s.read()
    
        with open(structure_file, "w") as o:
            o.write(line1)
            o.write(newline2)
            o.write(rest)    