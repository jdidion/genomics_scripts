#!/usr/bin/env python
# Summarize multiple FastQC reports in a single table, with links to the source reports.
# The input is a single directory, under which there may be from one to many library directories,
# and one to many FastQC reports in each library directory. If FastQC reports are zipped, they
# will be unzipped.
# Author: John Didion (jdidion@email.unc.edu)

from argparse import ArgumentParser
from csv import reader
from glob import glob
from mako.template import Template
import os
import sys

ICONS = dict(PASS="tick.png", WARN="warning.png", FAIL="error.png")

def read_summary_file(path):
    with open(path, "rU") as f:
        return tuple(row[0] for row in reader(f, delimiter="\t"))

class Library(object):
    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.reports = dict()
    
    def __cmp__(self, other):
        return cmp(self.name, other.name)
        
    def add_report_dir(self, name, path):
        self.reports[name] = Report(name, path)
    
    def add_report_zip(self, name, lib_path, zip_path):
        from zipfile import ZipFile
        with ZipFile(zip_path) as z:
            z.extractall(lib_path)
        self.add_report_dir(name, os.path.join(lib_path, name))
    
    def has_report(self, name):
        return name in self.reports
    
    def iter_reports(self):
        for name in sorted(self.reports.keys()):
            yield (name, self.reports[name])
        
class Report(object):
    def __init__(self, name, path):
        self.name = name
        self.display_name = name[0:(len(name)-7)] # strip _fastqc from name
        self.path = path
        self.summary = None
    
    def get_summary(self):
        if self.summary is None:
            summary_file = os.path.join(self.path, "summary.txt")
            self.summary = read_summary_file(summary_file)
        return self.summary
    
    def iter_summary(self):
        summary = self.get_summary()
        index = self.get_index_path()
        for i in xrange(0, len(summary)):
            icon = os.path.join(self.path, "icons", ICONS[summary[i]])
            link = "{0}#M{1}".format(index, i)
            yield (icon, link)
        
    def get_index_path(self):
        return os.path.join(self.path, "fastqc_report.html")

def main():
    parser = ArgumentParser()
    parser.add_argument("-o", "--outfile", default=None, 
        help="Alternative path to output summary file.")
    parser.add_argument("-t", "--template_file", default="fastqc_summary.html",
        help="Path to template file; defaults to the current working directory.")
    parser.add_argument("--nounzip", action="store_true", default=False,
        help="Do not unzip any reports.")
    parser.add_argument("dirs", nargs="+", help="Directory(ies) containing reports.")
    args = parser.parse_args()
    
    libraries = set()
    
    for lib_path in (p for d in args.dirs for p in glob(d)):
        if os.path.isdir(lib_path):
            lib_entry = os.path.basename(lib_path)
            lib = Library(lib_entry, lib_path)
            libraries.add(lib)
            
            for rep_entry in os.listdir(lib_path):
                rep_path = os.path.join(lib_path, rep_entry)
                
                if os.path.isdir(rep_path) and not lib.has_report(rep_entry):
                    lib.add_report_dir(rep_entry, rep_path)
                elif os.path.isfile(rep_path) and not args.nounzip:
                    rep_name, ext = os.path.splitext(rep_entry)
                    if ext == ".zip" and not lib.has_report(rep_name):
                        lib.add_report_zip(rep_name, lib_path, rep_path)
    
    outfile = args.outfile or os.path.join(args.dir, "summary.html")
    with open(outfile, 'w') as f: 
        f.write(Template(filename=args.template_file).render(dir=args.dir, libs=sorted(libraries)))
    
if __name__ == "__main__":
    main()
