import bamnostic as bs
from pathlib import Path
import numpy as np

from numba import jit, njit

home = Path(__file__).parent
decode_map = np.array(['A', 'C', 'G', 'T', 'N'])

class utils:

    def io_monster(mode, n_rows=0):

        def read_encode_bam():
            seqs = [read.seq for read in bs.AlignmentFile(bam_path)]
            lengths = np.array([len(seq) for seq in seqs], dtype=np.int32)                                 #store the indices of the reads so we can retrieve them later
            all_bytes = np.frombuffer("".join(seqs).encode(), dtype=np.uint8)   #flatten the whole genome into a byte array (strings -> bytes)
            encode_map = utils.init_base_encoding_map()
            seqs = utils.encode_sequences(all_bytes, encode_map)                 #use the fast numba function to encode the bytes as nucleotides (bytes -> encodings)
            indptr = np.zeros(len(lengths) + 1, dtype=np.int32)             #build indptr (map of indices in flat array)
            indptr[1:] = np.cumsum(lengths)                                 #get the fenceposts of all sequences in the flat array
            np.savez_compressed(npz_path, seqs=seqs, indptr=indptr)


        npz_name = "processed_data.npz"
        bam_path = home / "data" / "SRR9090854.subsampled_5pct.bam"
        npz_path = home / "data" / npz_name

        if mode == "test":                                                  #use a test file if testmode specified
            print("loading bam test data query sequences")
            bam = bs.AlignmentFile(bs.example_bam, 'rb')
            seqs = [seq.query_sequence for seq in bam]
        else:                                                               #If n_rows != 0, use a text file containing the specified number of rows
            if not npz_path.exists():                                       #text file does not exist
                print(f"{npz_name} does not exist, get ready for bammage. reading from bam.")
                read_encode_bam()
            else:
                print(f"{npz_name} exists, reading {n_rows} rows from path")
            npz = np.load(npz_path)
            print("data ready :)")
            return npz["seqs"], npz["indptr"]
        return seqs, None
    
    def init_base_encoding_map():
        base_map = np.zeros(256, dtype=np.int8)
        base_map[ord('A')] = 0
        base_map[ord('C')] = 1
        base_map[ord('G')] = 2
        base_map[ord('T')] = 3
        base_map[ord('N')] = 4
        return base_map
    
    def decode_sequence(seq):
        return "".join(decode_map[seq])

    
    def seq_to_array(seq, base_map):
        seq_array = np.frombuffer(seq.encode(), dtype=np.int8)
        indices = base_map[seq_array]
        return indices
    
    @njit
    def encode_sequences(all_bytes, base_map):
        mdata = np.empty(len(all_bytes), dtype=np.int8)
        for i in range(len(all_bytes)):
            mdata[i] = base_map[all_bytes[i]]
        return mdata
    
    @njit
    def fast_init(flat_seqs, n_rows, indptr, k):
        motifs = np.zeros((n_rows, k), dtype=np.int8)                   #init array for motif storage
        idxs = np.zeros(n_rows, dtype=np.int32)                         #init array for motif indices
        for i in range(n_rows):
            idx = np.random.randint(indptr[i], indptr[i+1] - k + 1)     #choose random index
            motifs[i] = flat_seqs[idx:idx+k]                            #get the motif
            idxs[i] = idx + indptr[i]                                   #save the chosen motif index at its position within that sequence
        return motifs, idxs
    
    @njit
    def build_pfm_fast(motifs, k, n_rows, n_bases=5):
        pfm = np.zeros((k, n_bases), dtype=np.int32)
        for i in range(n_rows):
            for j in range(k):
                pfm[j, motifs[i, j]] += 1
        return pfm


        