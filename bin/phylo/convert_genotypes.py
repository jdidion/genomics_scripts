#!/usr/bin/env python
"""
Create an interleaved Phylip or Nexus input file given a genotype file in the format::

    SNP_ID, REF, VAR, strain1, strain2...
    JAX123, A, C, 1, 3...
    
Genotypes can be in one of several formats (specified by the -g option):
* nuc: A/C/G/T/H/N
* allele: A=REF, B=VAR, H, N
* num: 1=REF, 2=H, 3=VAR (anything else is N)
* es (EIGENSTRAT format): 0=VAR, 1=HET, 2=REF, 9=N
* dip: diploid genotypes (AA/CC/GG/TT/NN/AC...)
* mac: minor allele count: -1=N, 0=hom for major allele, 1=het, 2=hom for minor allele

The genotype file can have additional columns, but the --ref_column, --var_column and 
--first_strain_column options must be specified.

If the number of SNPs in the source file is known, it should be supplied with the -m argument, 
otherwise the program will have to scan the entire file to determine the number of markers before 
reading it again to do the conversion.
"""
from collections import OrderedDict
import csv
from functools import partial
from math import ceil, log10
import os
import re
import sys
from cStringIO import StringIO

import jpd.util.cl
from jpd.util.io import write_dict
from jpd.util.gnet import nt

def main():
    def add_args(parser):
        parser.add_argument("-d", "--delimiter", metavar="CHAR", default=",", help="Input file delimiter")
        parser.add_argument("-f", "--output_format", metavar="FORMAT", choices=["phylip","relaxed","nexus","fasttree"],
            default="phylip", help="Output format.")
        parser.add_argument("-g", "--genotype_format", metavar="FORMAT", 
            choices=["nuc", "allele", "num", "es", "dip", "mac"], default="num", help="Genotype format.")
        parser.add_argument("-m", "--markers", action="extend", type="int_list", metavar="N",
            help="Comma-separated list of number of markers in each input file "\
                "(calculated if not supplied).")
        parser.add_argument("-n", "--sample_name_length", type=int, metavar="CHARS", default=10,
            help="Max number of characters in sample name.")
        parser.add_argument("-s", "--samples", action="extend", type="str_list", metavar="INDICES",
            help="List of strain names or indices to include. All strains are included if omitted.")
        parser.add_argument("--ref_column", type=int, metavar="NUM", default=1,
            help="Column of reference base.")
        parser.add_argument("--var_column", type=int, metavar="NUM", default=2,
            help="Column of variant base.")
        parser.add_argument("--first_strain_column", type=int, metavar="NUM", default=3,
            help="First column with genotype data.")
        parser.add_argument("--map_name_to_num", type='writeable_file', metavar="FILE", default=None,
            help="For phylip files, replace sample names with numbers and write a file with number=name.")
        parser.add_argument("--infer_alleles", action="store_true", default=False,
            help="Guess the reference and variant alleles; genotype_format must be 'nuc'.")
        parser.add_argument("--sequential", action="store_true", default=False,
            help="Output data in sequential rather than interleaved format.")
        parser.add_argument("--classes", default=None,
            help="If specified, output will be in phylip discrete character format. Classes are "\
                 "separated by |. Each class can have one or more genotype value, separated by commas.")
        parser.add_argument("--discrete_chars", type='str_list', default=['?','0','1','2','3','4','5','6','7'],
            help="Characters used to represent discrete classes. Only used if the classes parameter "\
                 "is specified.")
        parser.add_argument("infiles", action="extend", type="readable_file", metavar="FILE", 
            nargs="*", help="Input genotype file(s).")
        parser.add_argument("outfile", type="writeable_file", metavar="FILE", help="Output file.")
    
    args = util.cl.parse(add_args)
    
    assert not args.infer_alleles or args.genotype_format == 'nuc'
    
    ninput = len(args.infiles)
    total_markers = args.markers or []
    assert len(total_markers) == 0 or len(total_markers) == ninput

    if not total_markers:
        for f in args.infiles:
            with open(f) as infile:
                for j, l in enumerate(infile): pass
            total_markers.append(j)

    name_len = args.sample_name_length
    pad = 1
    if args.output_format == 'phylip':
        name_len = 10
        pad = 0
    
    linelen = max(80, name_len+pad) # TODO make configurable
    if args.genotype_format == 'dip' and (linelen - name_len - pad) % 2 > 0:
        linelen += 1
    
    classes = None
    if args.classes is not None:
        ndisc = len(args.discrete_chars)
        classes = {}
        # work around limitation of argparse that arguments that don't look like numbers
        # can't begin with a dash
        if args.classes[0] == '\\':
            args.classes = args.classes[1:]
        for i,cl in enumerate(args.classes.split('|')):
            if i >= ndisc:
                sys.exit("Phylip supports a maximum of %d discrete classes" % ndisc)
            for ch in cl.split(','):
                classes[ch] = args.discrete_chars[i]
    
    out = open(args.outfile, 'w')
    curpos = None
    header = None
    names = None
    samples = None
    seqs = []
    
    delimiter = args.delimiter
    if delimiter == "TAB":
        delimiter = "\t"
    
    for file_idx, f in enumerate(args.infiles):
        with open(f, 'rU') as infile:
            reader = csv.reader(infile, delimiter=delimiter)
            fheader = reader.next()
            markers = 0
            
            if header is None:
                header = fheader[args.first_strain_column:]
                
                if args.samples:
                    samples = []
                    for s in args.samples:
                        if re.match('^\d+$', s):
                            samples.append(int(s))
                        elif re.match('^\d+-\d+$', s):
                            rng = map(int, s.split("-"))
                            samples.extend(xrange(rng[0]-1, rng[1]))
                        elif s in header:
                            samples.append(header.index(s))
                        else:
                            raise Exception("Invalid sample name %s" % s)
                    header = [header[i] for i in samples]
                
                nmarkers = sum(total_markers)
                if args.genotype_format == "dip":
                    nmarkers *= 2
                hdata = (len(header), nmarkers)
                if args.output_format in ('phylip','relaxed','fasttree'):
                    print >>out, "%i %i" % hdata
                elif args.output_format == 'nexus':
                    print >>out, "begin data;"
                    print >>out, "dimensions ntax=%i nchar=%i;" % hdata
                    print >>out, "format datatype=dna interleave=yes gap=-;"
                    print >>out, "matrix"
                
                # generate sample names
                if args.map_name_to_num:
                    zpad = int(ceil(log10(len(header)))) + 1
                    if name_len is None:
                        name_len = zpad
                    names = [str(i).rjust(zpad, '0') for i in xrange(1,len(header)+1)]
                    write_dict(args.map_name_to_num, OrderedDict(zip(names, header)))
                else:
                    if name_len is None:
                        name_len = max(len(s) for s in header)
                    names = [re.sub('\W', '_', s[0:name_len]) for s in header]
                names = [n.ljust(name_len+pad,' ') for n in names]

                # initialize buffers
                for i,n in enumerate(names):
                    seqs.append(StringIO())
                    seqs[i].write(n)
            else:
                assert header == fheader[args.first_strain_column:], "Headers do not match"
            
            if curpos is None:
                curpos = name_len + pad
            
            for row in reader:
                if markers == total_markers[file_idx]:
                    break
                    
                curpos += 2 if args.genotype_format == 'dip' else 1
                markers += 1
                
                row_genos = row[args.first_strain_column:]
                
                genos = None
                if classes is not None:
                    genos = classes
                elif args.genotype_format != 'dip':
                    if args.infer_alleles:
                        alleles = set(row_genos) - set(('H','N'))
                        if len(alleles) == 0:
                            # This is a hack to deal with SNPs with all H/N genotypes
                            ref = 'C'
                            var = 'G'
                        elif len(alleles) == 1:
                            ref = var = alleles.pop()
                        elif len(alleles) == 2:
                            ref, var = alleles
                        else:
                            raise Exception("Cannot infer alleles for row %s" % row[0])
                    else:
                        ref = row[args.ref_column]
                        var = row[args.var_column]
                        
                    het = nt.get_or(ref, var) if ref and var else ''
                    if args.genotype_format == 'num':
                        genos = { '1': ref, '2': het, '3': var }
                    elif args.genotype_format == 'es':
                        genos = { '2': ref, '1': het, '0': var }
                    elif args.genotype_format == 'allele':
                        genos = { 'A': ref, 'H': het, 'B': var }
                    elif args.genotype_format == 'nuc':
                        genos = { ref: ref, 'H': het, var: var }
                
                if samples:
                    row_genos = [row_genos[i] for i in samples]
                for i,g in enumerate(row_genos):
                    seqs[i].write(genos.get(g, 'N') if genos else g)

                # restrict lines to linelen characters
                if not args.sequential and curpos % linelen == 0:
                    if args.output_format in ('nexus','fasttree'):
                        curpos += name_len + pad
                    
                    for i,s in enumerate(seqs):
                        print >>out, s.getvalue()
                        seqs[i] = StringIO()
                        if args.output_format in ('nexus','fasttree'):
                            seqs[i].write(names[i])
                    
                    # a newline must separate sequence blocks
                    print >>out 
    
    if args.sequential or curpos % linelen > 0:
        for s in seqs:
            print >>out, s.getvalue()
    
    if args.output_format == 'nexus':
        print >>out, ";\nend;"
        
    out.close()

if __name__ == '__main__':
    main()
