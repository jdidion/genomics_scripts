#!/usr/bin/env python

"""
Convert plink data to sequence formats (phylip, nexus, fasta).
"""

from jpd.gnet.plink import PlinkReader
from jpd.gnet.phylo import *
from jpd.gnet import nt
from jpd.util.cl import parse

GENOTYPE_FORMAT = dict(
    haploid=lambda A, B, N: { 0 : A, 1 : N, 2 : B, 3 : N },
    IUPAC=lambda A, B, N: { 0 : A, 1 : nt.get_or(A, B), 2 : B, 3 : N },
    diploid=lambda A, B, N: {
        0 : "{0}{0}".format(A),
        1 : "{0}{1}".format(*sorted((A, B))),
        2 : "{0}{0}".format(B),
        3 : "{0}{0}".format(N)
    }
)

def main():
    def add_args(parser):
        parser.add_argument("-b", "--bfile", metavar="PATH", required=True,
            help="PLINK binary file prefix.")
        parser.add_argument("-c", "--chrm", metavar="CHRM", default=None,
            help="Chromosome to process (or whole genome if None).")
        parser.add_argument("-f", "--families", action="append", default=None,
            help="Family IDs to keep (or all families if None).")
        parser.add_argument("-g", "--genotype-format", metavar="FORMAT",
            choices=("haploid", "diploid", "IUPAC"), default="haploid",
            help="Genotype format; haploid = A/C/G/T/N; diploid = output two characters for "\
                 "each site; IUPAC = output a single character using IUPAC codes for "\
                 "ambiguous sites.")
        parser.add_argument("-m", "--map-name-to-num", action="store_true", default=False,
            help="For phylip files, replace sample names with numbers and write a "\
                 "file with number=name.")
        parser.add_argument("-o", "--outfile", type="writeable_file", metavar="FILE", 
            help="Output file prefix (full path without file extension).")            
        parser.add_argument("-O", "--output-format", metavar="FORMAT", 
            choices=("phylip", "relaxed", "nexus", "fasttree", "fasta"), default="phylip", 
            help="Output format.")
        parser.add_argument("-s", "--samples", action="extend", type="str_list", metavar="INDICES",
            help="List of strain names or indices to include. All strains are included if omitted.")
        parser.add_argument("--sequential", action="store_true", default=False,
            help="Output data in sequential rather than interleaved format.")
    
    args = parse(add_args)
    
    reader = PlinkReader(args.bfile, args.chrm, args.families, args.samples)
    sample_ids = reader.sample_ids
    outfile = args.outfile or args.bfile
    
    nchar = None
    if args.map_name_to_num :
        nchar = len(str(len(sample_ids)))
        new_sample_ids = list(str(i).zfill(nchar) for i in xrange(len(sample_ids)))
        map_file = "{0}.txt".format(outfile)
        from csv import writer
        with open(map_file, "w") as o:
            w = writer(o)
            for row in zip(new_sample_ids, sample_ids):
                w.writerow(row)
        sample_ids = new_sample_ids
    
    geno_map = GENOTYPE_FORMAT[args.genotype_format]
    
    if args.output_format == "fasta":
        writer = FastaWriter("{0}.fasta".format(outfile), sample_ids, geno_map)
    
    else:
        num_loci = reader.num_loci
        if args.genotype_format == "diploid":
            num_loci *= 2

        if args.output_format in ("phylip", "relaxed"):
            file_format = PhylipFormat(sample_ids, num_loci, 
                strict=args.output_format == "phylip", 
                interleaved=not args.sequential)
                
        elif args.output_format == "nexus":
            file_format = NexusFormat(sample_ids, num_loci, 
                interleaved=not args.sequential)
        
        elif args.output_format == "fasttree":
            file_format = FastTreeFormat(sample_ids, num_loci,
                interleaved=not args.sequential)
        
        writer = PhyloWriter(file_format, "{0}.phy".format(outfile), geno_map)
    
    try:
        for locus, genotypes in reader:
            writer.add_snp(locus, genotypes)

    finally:
        reader.close()
        writer.finish()

if __name__ == "__main__":
    main()
