import bamnostic as bs
from pathlib import Path
import numpy as np

from numba import jit, njit

home = Path(__file__).parent
decode_map = np.array(['A', 'C', 'G', 'T', 'N'])

class utils:

    def io_monster(mode, n_rows=0):

        def read_bam_subset():
            seqs = []
            data = bs.AlignmentFile(bam_path)
            i = 0
            for read in data:
                seqs.append(read.seq)
                if i == n_rows:
                    break
                i += 1
            return seqs
        
        def write_txt_file(seqs):
            with open(txt_path, "w", encoding="utf-8") as f:        
                for line in seqs:
                    f.write(f"{line}\n")
        
        def read_txt_file():
            seqs = []
            with open(txt_path, "r", encoding="utf-8") as f:
                for line in f.readlines():
                    seqs.append(line.strip())
            return seqs

        bam_path = home / "data" / "SRR9090854.subsampled_5pct.bam"
        txt_path = home / "data" / "converted_data.txt"

        if mode == "test":                                                  #use a test file if testmode specified
            print("loading bam test data query sequences")
            bam = bs.AlignmentFile(bs.example_bam, 'rb')
            seqs = [seq.query_sequence for seq in bam]
        elif n_rows == 0:                                                   #READ THE WHOLE THING
            seqs = [read.seq for read in bs.AlignmentFile(bam_path)]
        else:                                                               #If n_rows != 0, use a text file containing the specified number of rows
            if not txt_path.exists():                                       #text file does not exist
                print(f"converted_data.txt does not exist, get ready for bammage. reading {n_rows} rows from bam.")
                seqs = read_bam_subset()
                write_txt_file(seqs)
            else:
                print(f"converted_data.txt exists, reading {n_rows} rows from path")
                seqs = read_txt_file()
                if len(seqs) != n_rows:
                    print(f"row mismatch; reading from bam. CANCEL IF THE n_rows IS WRONG")
                    seqs = read_bam_subset()
                    write_txt_file(seqs)

        print("data ready :)")
        return seqs
    
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
    def fast_init(flat_seqs, n_rows, lengths, k):
        motifs = np.zeros((n_rows, k), dtype=np.int8)                   #init array for motif storage
        indptr = np.zeros(len(lengths) + 1, dtype=np.int32)             #build indptr (map of indices in flat array)
        indptr[1:] = np.cumsum(lengths)                                 #get the fenceposts of all sequences in the flat array
        idxs = np.zeros(n_rows, dtype=np.int32)                         #init array for motif indices
        for i in range(n_rows):
            idx = np.random.randint(indptr[i], indptr[i+1] - k + 1)     #choose random index
            motifs[i] = flat_seqs[idx:idx+k]                            #get the motif
            idxs[i] = idx + indptr[i]                                   #save the chosen motif index at its position within that sequence
        return motifs, idxs, indptr
    
    @njit
    def build_pfm_fast(motifs, k, n_rows, n_bases=5):
        pfm = np.zeros((k, n_bases), dtype=np.int32)
        for i in range(n_rows):
            for j in range(k):
                pfm[j, motifs[i, j]] += 1
        return pfm


        