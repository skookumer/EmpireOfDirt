import numpy as np
from numba import njit
from utils import utils



class sequence_box:

    def __init__(self, indptr: np.array, seq_array: np.array, motif_indices: np.array, motifs: np.array, k: int, n_bases=5):
        self.indptr = indptr
        self.seqs = seq_array
        self.n_rows = len(indptr) - 1                               #subtract 1 because indptr is +1 longer than the array (the last fence post)
        self.midx = motif_indices
        self.motifs = motifs
        self.k = k
        self.n_bases = n_bases
    
    def __len__(self):
        return self.n_rows
    
    def __getitem__(self, i):
        if i < 0 or i >= self.n_rows:
            raise IndexError("INDEX OUT OF RANGE :(")
        x = self.indptr[i]
        y = self.indptr[i+1]
        return self.seqs[x:y]

    def __iter__(self):
        for i in range(self.n_rows):
            yield self.__getitem__(i)

    def init_bg(self):
        total_freqs = np.bincount(self.seqs)
        pfm = utils.build_pfm_fast(self.motifs, self.k, self.n_rows, self.n_bases)
        motif_freqs = pfm.sum(axis=0)
        print(total_freqs)
        print(motif_freqs)
        if total_freqs.shape[0] != motif_freqs.shape[0]:
            raise IndexError("array mismatch")
        return total_freqs - motif_freqs
    
    def init_pfm(self):
        return utils.build_pfm_fast(self.motifs, self.k, self.n_rows, self.n_bases)
