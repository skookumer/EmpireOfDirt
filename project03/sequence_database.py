import numpy as np
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

    def get_bg(self):
        print("getting background frequencies\n")
        total_freqs = np.bincount(self.seqs)
        pfm = utils.build_pfm_fast(self.motifs, self.k, self.n_rows, self.n_bases, slice_last_row=False)
        motif_freqs = pfm.sum(axis=1)
        if total_freqs.shape[0] != motif_freqs.shape[0]:
            raise IndexError("array mismatch :((")
        return total_freqs - motif_freqs
    
    def get_pfm(self, to_mask=None):
        print("building frequency matrix\n")
        if to_mask is not None:
            return utils.build_pfm_fast(self.motifs, self.k, self.n_rows, self.n_bases)
        mask = np.ones(shape=self.n_rows, dtype=np.bool_)
        mask[to_mask] = False
        return utils.build_pfm_fast(self.motifs[mask], self.k, self.n_rows - 1, self.n_bases)

    def get_str_list_format_motifs(self):
        return [utils.decode_sequence(entry) for entry in self.motifs]

    def get_str_list_format_seqs(self):
        output = []
        for i in range(len(self.indptr) - 1):
            x = self.indptr[i]
            y = self.indptr[i + 1]
            output.append(utils.decode_sequence(self.seqs[x:y]))
        return output