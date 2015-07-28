#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# binbam: Partition a BAM file by template length.
#

from __future__ import print_function
import argparse
import getpass
import os
import sys

import pysam


def parse_intervals(intervals_string):
    intervals = [interval_string.split('-') for interval_string in intervals_string.split(',')]

    cleaned_intervals = []
    for interval in intervals:
        if len(interval) != 2:
            raise argparse.ArgumentTypeError("Intervals are specified with exactly two integers.")

        start, end = [int(s) for s in interval]

        if not start:
            start = 1
            print("Interval %s start omitted; using 1." % '-'.join(interval), file=sys.stderr)
        elif start < 1:
            raise argparse.ArgumentTypeError("Intervals must consist of positive integers.")
        else:
            start = int(start)

        if not end:
            end = float('inf')
            print("Interval %s end omitted; using INFINITY!" % '-'.join(interval), file=sys.stderr)
        elif end < 1:
            raise argparse.ArgumentTypeError("Intervals must consist of positive integers.")
        else:
            end = int(end)

        if start >= end:
            new_interval = (end, start)
            print("Interval %s specified backward; correcting to %s" % ('-'.join(interval), '-'.join(str(n) for n in new_interval)), file=sys.stderr)
        else:
            new_interval = (start, end)

        cleaned_intervals.append(new_interval)

    cleaned_intervals = sorted(cleaned_intervals)
    last_start, last_end = 0, 0
    for interval in cleaned_intervals:
        start, end = interval
        if start <= last_end:
            raise argparse.ArgumentTypeError("Interval %d-%d overlaps %d-%d." % (start, end, last_start, last_end))
        last_start, last_end = start, end

    return cleaned_intervals


def select_bin(length, intervals):
    fits = [start <= abs(length) <= end for start, end in intervals]
    try:
        return fits.index(True)
    except ValueError:
        return None


def make_bin_name(source_filename, interval):
    base = os.path.basename(source_filename)
    return base.replace('.bam', '.bin.%d.%d.bam' % interval)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Partition a BAM file by template length, creating a bin BAM file for each interval given.')
    parser.add_argument('-i', '--intervals', dest='intervals', type=parse_intervals, required=True, help='A comma-separated list of intervals defined as two positive integers separated by a hyphen, e.g. "25-75,125-175". It is an error for them to overlap.')
    parser.add_argument('source', help='The BAM file to partition.')
    args = parser.parse_args()

    source = pysam.AlignmentFile(args.source, "rb")

    # args.intervals is a sorted list of intervals as tuples. make a
    # file for each interval and keep them in a list in the same order
    bin_files = []
    for interval in args.intervals:
        bin_files.append(pysam.AlignmentFile(make_bin_name(args.source, interval), "wb", header=source.header))

    for i, record in enumerate(source):
        length = record.template_length
        bin = select_bin(length, args.intervals)
        if bin is not None:
            bin_files[bin].write(record)

    for f in bin_files:
        f.close()
