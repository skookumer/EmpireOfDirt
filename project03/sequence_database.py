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
        self.selected_motif = None
        self.frozen = np.zeros(shape=self.n_rows, dtype=bool)
        self.sampling_pool = np.zeros(shape=self.n_rows, dtype=bool)
    
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
    
    def get_pfm(self, to_mask=False):
        if not to_mask:
            return utils.build_pfm_fast(self.motifs, self.k, self.n_rows, self.n_bases)
        mask = self.frozen
        mask[self.selected_motif] = False
        motif_array = self.motifs[mask]
        return utils.build_pfm_fast(motif_array, self.k, len(motif_array), self.n_bases)

    def get_str_list_format_motifs(self):
        return [utils.decode_sequence(entry) for entry in self.motifs]

    def get_str_list_format_seqs(self):
        output = []
        for i in range(len(self.indptr) - 1):
            x = self.indptr[i]
            y = self.indptr[i + 1]
            output.append(utils.decode_sequence(self.seqs[x:y]))
        return output
    
    def select_random_motif(self):
        self.selected_motif = np.random.randint(0, len(self.motifs))
        return self.__getitem__(self.selected_motif)
    
    def update_motifs(self, motif):
        self.motifs[self.selected_motif] = motif

    def unfreeze_random(self, subsample_size):
        frozen_indices = np.where(self.frozen == 0)[0]
        if len(frozen_indices) > subsample_size:
            samp_size = subsample_size
        else:
            samp_size = len(frozen_indices)

        indices = np.random.choice(frozen_indices, size=samp_size, replace=False)
        self.frozen[indices] = True
        self.sampling_pool[indices] = True
    
    def reset_sampling_pool(self):
        self.sampling_pool = np.zeros(shape=self.n_rows, dtype=bool)

    def select_random_sampling_motif(self):
        sampling_indices = np.where(self.sampling_pool == 1)[0]
        self.selected_motif = np.random.choice(sampling_indices)
        return self.__getitem__(self.selected_motif)
    
    def split_subsamples(self, subsample_size):
        i = 0
        j = 0
        subsamples = []
        while i < self.n_rows:
            j = min(i + subsample_size, self.n_rows)
            offset = self.indptr[i]
            indptr = [value - offset for value in self.indptr[i: j]]
            motifs = self.motifs[i: j]
            motif_indices = self.midx[i: j]
            seq_array = self.seqs[self.indptr[i]: self.indptr[j]]
            subsamples.append(sequence_box(indptr, seq_array, motif_indices, motifs, self.k, self.n_bases))
            i += subsample_size
        return subsamples
    
    def check_remaining_frozen(self):
        if len(np.where(self.frozen == 0)[0]) > 5:
            return True
        return False
