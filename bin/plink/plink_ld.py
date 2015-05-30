#!/usr/bin/env python
"""
Command line interface to the plink LD functions in plink.py. The general workflow is:
1. partition_snps: Split a plink .bim file into "ld-snp-list" files (a list of SNP ids corresponding 
   to all the SNPs in a single window).
2. run plink: This uses the plink arguments --ld-window-r2 0 --inter-chr to perform all pairwise
   comparisons. Unless the --unfiltered argument is specified, it is assumed that the plink files
   have been filtered for uncalled and monomorphic genotypes.
3. process_ld_files: Parse the .ld files output from plink and convert them into component files
   that will be merged together in the next step.
4. merge_ld_files: Merge component files into summary files suitable for plotting and other analysis.
5. plot: Create a triangle heatmap plot.
"""

import fileinput
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.environ['SCI_HOME'], "common", "lib"))
import geno.plink
from geno.plink import *
from viz.ld import plot_ld_heatmap
from util.cl import parse, readable_file_group, delimited_macro
from util.fork import *
from util.io import mkdir
from util.misc import bash

FORK_COMMANDS = ["plink","split","blocks","block_ld"]
SERIAL_COMMANDS = ["snps","process","plot","local","process_local","prune_blocks", "tag_snps"]
OTHER_COMMANDS = ["make_blocks"]
COMMANDS = FORK_COMMANDS + SERIAL_COMMANDS + OTHER_COMMANDS

def call_make_blocks(args):
    make_blocks(*args)

def main():
    def add_args(parser):
        parser.add_argument("-b", "--bfile", type=readable_file_group(("bed", "bim", "fam")),
            default=None, help="Prefix of plink bfiles (bed, bim and fam).")
        parser.add_argument("-c", "--chromosomes", type=delimited_macro("chrm"), default="F",
            help="Set of chromosomes on which to execute the command(s)")
        parser.add_argument("-f", "--fork_mode", choices=("test","serial","thread","lsf"),
            default="thread", help="How to distribute jobs.")
        parser.add_argument("-i", "--interval_size", type=int, default=1000, metavar="bp",
            help="Bin size for interval statistics (in bp)")
        parser.add_argument("-k", "--keep_file", type="readable_file", metavar="FILE", default=None,
            help="Plink keep file (list of samples to include).")
        parser.add_argument("-s", "--plot_stats", action="extend_overwrite", type="str_list", 
            metavar="LIST", default=("mean","max","pct"),
            help="List of stats for which to create (mean, max, pct)")
        parser.add_argument("-S", "--summary_stat", action="extend_overwrite", 
            choices=("all", "unlinked"), default="unlinked",
            help="How to compute bin summary stats (use all markers or only unlinked markers)")
        parser.add_argument("-p", "--percentile", type=int, default=95, metavar="PCT",
            help="Percentile for r-squared statistics.")
        parser.add_argument("-r", "--r2_bin_size", type=float, default=0.01,
            help="Bin size for r-squared histogram.")
        parser.add_argument("-w", "--window_size", type=int, default=500, metavar="Kb",
            help="Window size (in kb)")
        parser.add_argument("--fork_opts", action="extend_dict", type="mapping_list", default={},
            help="Options specific to the fork mode.")
        parser.add_argument("--unfiltered", action="store_true", default=False,
            help="Assume data files are unfiltered and apply uncalled and maf filters.")
        parser.add_argument("--no_unlinked", action="store_true", default=False,
            help="Don't compute stats for unlinked markers.")
        parser.add_argument("--per_chromosome", action="store_true", default=False,
            help="Whether the data file has been split into one per chromosome.")
        parser.add_argument("--block_r2_thresholds", type=float, nargs=2, default=(0.2, 0.95),
            help="Lower and upper r-square thresholds for identifying haplotype blocks with 'make_blocks' command.")
        parser.add_argument("--output_format", choices=("ped","bed","tped"), default="ped",
            help="Output format for commands that produce PLINK data files.")
        parser.add_argument("--new_plink", default="/Users/johndidion/software/plink_mac/plink")
        parser.add_argument("--old_plink", default="plink")
        parser.add_argument("outdir", type="writeable_dir")
        parser.add_argument("commands", action="extend", nargs="+", choices=COMMANDS,
            help="Commands to run. If none are specified, all will be run.")
    
    ns = parse(add_args)
    
    geno.plink.PLINK_CMD = ns.new_plink
    geno.plink.PLINK_OLD = ns.old_plink
    
    if ns.bfile:
        bedfile, bimfile, famfile = ns.bfile
        bfile = os.path.splitext(bedfile)[0]
        fname = os.path.basename(bfile)

    window_size = ns.window_size * 1000
    commands = ns.commands if ns.commands else COMMANDS
    
    mkdir(ns.outdir, overwrite=False)
    
    if any(c in FORK_COMMANDS for c in commands):
        executor = get_executor(ns.fork_mode, ns.fork_opts)
    
    # partition the genome into windows and generate a list of SNPs in each window
    if "snps" in commands:
        # this command is not forked but is not time-consuming
        chr_file, win_file = partition_snps(bimfile, ns.outdir, window_size)
    else:
        win_file = os.path.join(ns.outdir, "windows.csv")
    
    # execute the plink --r2 command over snp windows (requires 'snps' command)
    if "plink" in commands:
        cmd_iter = pairwise_ld_command_iter(bfile, ns.outdir, win_file, apply_filters=ns.unfiltered)
        exec_shell(cmd_iter, executor, error_handler=reraise_error)
    
    # process the results of the 'plink' command
    if "process" in commands:
        # this command is not forked. it should be submitted to lsf if run on kure.
        process_ld_files(ns.outdir, ns.percentile, ns.interval_size, ns.r2_bin_size, unlinked=not ns.no_unlinked)
    
    # generate heatmap plots from the results of the 'process' command
    if "plot" in commands:
        bins = os.path.join(ns.outdir, "window_summary.csv")
        for s in ns.plot_stats:
            summary_stat = "{0}_{1}".format(ns.summary_stat, s)
            plot_ld_heatmap(
                os.path.join(ns.outdir, "{0}_r2_matrix.csv".format(s)), bins, 
                os.path.join(ns.outdir, "{0}_heatmap.pdf".format(s)), 
                window_size=window_size, summary_stat=summary_stat, tics=True)
    
    if "local" in commands:
        # TODO: execute commands for local LD
        pass
    
    if "process_local" in commands:
        process_local_ld_file(os.path.join(ns.outdir, "local.ld"), os.path.join(ns.outdir, "r2_hist.csv"))
    
    per_chrm = ns.per_chromosome
    
    # split a whole-genome data file into one file per chromosome
    if "split" in commands:
        cmd_iter = split(bfile, ns.outdir, chromosomes=ns.chromosomes, output_format=ns.output_format)
        exec_shell(cmd_iter, executor, error_handler=reraise_error)
        per_chrm = True
    
    # execute the plink --blocks command
    if "blocks" in commands:
        # this command is only forked if run on split files (1 per chr). it should be submitted to lsf if run on kure.
        if per_chrm:
            cmd_iter = per_chrm_iter(bfile, ns.outdir, format_ld_blocks_command, ns.chromosomes, apply_filters=ns.unfiltered)
            exec_shell(cmd_iter, executor, error_handler=reraise_error)
        else:
            cmd = format_ld_blocks_command(bfile, ns.outdir, apply_filters=ns.unfiltered)
            bash(cmd)

    if "block_ld" in commands:
        if per_chrm:
            path = os.path.join(ns.outdir, "**", "*.blocks.det")
        else:
            path = os.path.join(ns.outdir, "*.blocks.det")
        
        lines = fileinput.input(glob.glob(path))
        cmd_iter = block_ld(lines, bfile, ns.outdir)
        exec_shell(cmd_iter, executor, error_handler=reraise_error)
        
    if "prune_blocks" in commands:
        prune_blocks(ns.outdir, ns.chromosomes)
    
    if "make_blocks" in commands:
        def arg_iter():
            for chrm in ns.chromosomes:
                yield (chrm, os.path.join(ns.outdir, "chr{0}".format(chrm), fname))
        distribute(call_make_blocks, arg_iter(), mode=ns.fork_mode, **ns.fork_opts)
        
    if "tag_snps" in commands:
        bash(format_merge_tag_lists_command(ns.outdir, fname, ns.chromosomes))
    
if __name__ == "__main__":
    main()