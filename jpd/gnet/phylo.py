from cStringIO import StringIO
import re

class Buffer(object):
    def __init__(self, init=None, init_type="first"):
        self.buf = StringIO()
        self.length = 0
        self.init = None
        if init is not None:
            self.buf.write(init)
            self.length = len(init)
            if init_type == "all":
                self.init = init
    
    def append(self, s):
        self.buf.write(s)
        self.length += len(s)
    
    def get(self, nchars=None):
        val = self.buf.getvalue()
        val_len = self.length
        self.buf.close()
        self.buf = StringIO()
        self.length = 0
        if self.init is not None:
            self.buf.write(self.init)
            self.length = len(self.init)
        if nchars is not None and val_len > nchars:
            self.buf.write(val[nchars:])
            self.length += (val_len - nchars)
            val = val[0:nchars]
        return val

class PhyloWriter(object):
    def __init__(self, file_format, outfile, geno_map):
        self.file_format = file_format
        self.outfile = outfile
        self.fh = open(outfile, "w")
        self.geno_map = geno_map
        
        if file_format.header is not None:
            self.fh.write(file_format.header)
        
        self.seqs = list(
            Buffer(
                init=None if file_format.sample_ids is None else file_format.sample_ids[i], 
                init_type=file_format.write_samples
            )
            for i in xrange(file_format.num_samples)
        )
    
    @property
    def seq_len(self):
        return self.seqs[0].length
    
    def add_snp(self, locus, genotypes):
        geno_map = self.geno_map(locus.allele1, locus.allele2, self.file_format.missing)
        nchars = len(geno_map.values()[0])
        for i, g in enumerate(genotypes):
            g = geno_map[g]
            assert len(g) == nchars
            self.seqs[i].append(g)
        if self.file_format.interleaved and self.seq_len >= self.file_format.line_width:
            self.write_sequences()
        
    def write_sequences(self):
        if self.seq_len > 0:
            nchars = min(self.seq_len, self.file_format.line_width)
            for s in self.seqs:
                temp = s.get(nchars)
                assert len(temp) == nchars
                self.fh.write(temp)
                self.fh.write("\n")
            self.fh.write("\n")
    
    def finish(self):
        self.write_sequences()
        if self.file_format.footer is not None:
            self.fh.write(file_format.footer)
        self.close()
    
    def close(self):
        if self.fh is not None:
            self.fh.close()
            self.fh = None

class FileFormat(object):
    def __init__(self, num_samples):
        self.header = None
        self.footer = None
        self.missing = "N"
        self.interleaved = True
        self.line_width = 80
        self.num_samples = num_samples
        self.sample_ids = None
        self.sample_id_len = 0
        self.write_samples = "first"

class PhylipFormat(FileFormat):
    def __init__(self, sample_ids, num_loci, strict=True, interleaved=True):
        FileFormat.__init__(self, len(sample_ids))
        
        self.header = "{0} {1}\n".format(len(sample_ids), num_loci)
        self.interleaved = interleaved
        
        # sanitize sample IDs
        sample_ids = map(lambda i: re.sub('\W', '_', i), sample_ids)
    
        # add sample names to the first block
        if strict:
            self.sample_ids = list(i[0:10].ljust(10, ' ') for i in sample_ids)
        
        else:
            self.sample_id_len = max(len(i) for i in sample_ids) + 1
            self.sample_ids = list(i.ljust(self.sample_id_len) for i in sample_ids)

class NexusFormat(FileFormat):
    def __init__(self, sample_ids, num_loci, interleaved=True):
        FileFormat.__init__(self, len(sample_ids))
        self.header = "begin data;\n"\
            "dimensions ntax={0} nchar={1};\n"\
            "format datatype=dna interleave={2} gap=-;\n"\
            "matrix".format(len(sample_ids), num_loci, "yes" if interleaved else "no")
        self.footer = ";\nend;"
        self.interleaved = interleaved
        self.sample_id_len = max(len(i) for i in sample_ids) + 1
        self.sample_ids = list(i.ljust(self.sample_id_len) for i in sample_ids)
        self.write_samples = "all"

class FastTreeFormat(FileFormat):
    def __init__(self, sample_ids, num_loci, interleaved=True):
        FileFormat.__init__(self, len(sample_ids))
        self.header = "{0} {1}\n".format(len(sample_ids), num_loci)
        self.interleaved = interleaved
        self.sample_id_len = max(len(i) for i in sample_ids) + 1
        self.sample_ids = list(i.ljust(self.sample_id_len) for i in sample_ids)
        self.write_samples = "all"

class FastaWriter(object):
    def __init__(self, outfile, sample_ids, geno_map, line_width=80):
        self.outfile = outfile
        self.sample_ids = sample_ids
        self.geno_map = geno_map
        self.line_width = 0
        self.seqs = list(StringIO() for i in xrange(len(sample_ids)))
    
    def add_snp(self, locus, genotypes):
        geno_map = self.geno_map(locus.allele1, locus.allele2, "N")
        nchars = len(genotypes.values()[0])
        for i, g in enumerate(genotypes):
            g = self.geno_map[g]
            assert len(g) == nchars
            self.seqs[i].write(g)
        
    def finish(self):
        def wrap(s, width):
            nchars = len(s)
            i = 0
            while i < nchars:
                x = min(i+width, nchars)
                temp = s[i:x]
                i += x
                yield temp
        
        with open(outfile, "w") as fasta:
            for i in xrange(len(self.sample_ids)):
                fasta.write("> {0}\n".format(self.sample_ids[i]))
                fasta.writelines(wrap(self.seqs[i].getvalue(), self.line_width))

    def close(self):
        pass
