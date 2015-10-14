from plinkio import plinkfile
from itertools import izip
from jpd.util.collections import index_map

class PlinkReader(object):
    """
    Iterate over data from PLINK binary files (BED/BIM/FAM). Each row is returned as a 
    two-element tuple: ( pos, (geno1, geno2, ...) ).
    """
    def __init__(self, path, chrm=None, fids=None, iids=None):
        self.chrm = chrm
        self.handle = plinkfile.open(path)
        if not self.handle.one_locus_per_row():
            raise Exception("This script requires that SNPs are rows and samples columns.")
        
        samples = self.handle.get_samples()
        self._subset_idxs = None
        
        if fids is None and iids is None:
            self._samples = dict((s.iid, s) for s in samples)
        
        else:
            if fids is not None:
                fids = set(fids)
            if iids is not None:
                iids = set(iids)
            def keep(s):
                return (fids is None or s.fid in fids) and (iids is None or s.iid in iids)
            self._subset_idxs = set(i for i, sample in enumerate(samples) if keep(sample))
            self._samples = index_map(samples[i].iid for i in self._subset_idxs)
        
        self._loci = self.handle.get_loci()
        self._iter = izip(self._loci, self.handle)
    
    @property
    def num_samples(self):
        return len(self._samples)
    
    def has_sample(self, name):
        return name in self._samples
    
    @property
    def sample_ids(self):
        return self._samples.keys()
    
    def get_sample_index(self, name):
        return self._samples[sample]
    
    @property
    def num_loci(self):
        return len(self.handle.get_loci())
    
    def __iter__(self):
        return self
    
    def next(self):
        try:
            locus, genotypes = self._iter.next()
            while self.chrm is not None and str(locus.chromosome) != self.chrm:
                locus, genotypes = self._iter.next()
            if self._subset_idxs is not None:
                genotypes = list(g for i, g in enumerate(genotypes) if i in self._subset_idxs)
            return (locus, genotypes)

        except StopIteration:
            self.close()
            raise
    
    def close(self):
        self.handle.close()
