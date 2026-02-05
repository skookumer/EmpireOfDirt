

import random
import numpy as np
import bamnostic as bs
import seqlogo as sl
from pprint import pprint

#import function for building sequence motif & idenfitying seqs matching to motif
from data_readers import *
import seq_ops
import motif_ops

from joblib import parallel, delayed

from utils import utils
from sequence_database import sequence_box
from collections import Counter, defaultdict

bases = ["A", "T", "C", "G", "N"]
encode_map = utils.init_base_encoding_map()



def GibbsMotifFinder(seqs=None, k=6, n_rows=3083497, mode="norm", max_iter=1024, seed=None, toprint=True):
    '''
    Function to find a pfm from a list of strings using a Gibbs sampler
    
    Args: 
        seqs (str list): a list of sequences, not necessarily in same lengths
        k (int): the length of motif to find
        seed (int, default=None): seed for np.random

    Returns:
        pfm (numpy array): dimensions are 4xlength
    '''

    '''
    -----------------------------------------------INITIALIZATION---------------------------------------
    '''

    def initialize_motifs(seqs, k):                       #random pick of k-length sequences for each line of DNA as Motifs
        if toprint:
            print("randomly choosing motifs from each read")
        lengths = np.array([len(seq) for seq in seqs], dtype=np.int32)                                 #store the indices of the reads so we can retrieve them later
        all_bytes = np.frombuffer("".join(seqs).encode(), dtype=np.uint8)   #flatten the whole genome into a byte array (strings -> bytes)
        seqs = utils.encode_sequences(all_bytes, encode_map)                 #use the fast numba function to encode the bytes as nucleotides (bytes -> encodings)
        motifs, idxs, indptr = utils.fast_init(seqs, n_rows, lengths, k)             #use another fast numba function to randomly pick motifs and save the indices
        seq_box = sequence_box(indptr, seqs, idxs, motifs, k)                                                  #store the flat motif array with the indices in a python object for easy retrieval
        if toprint:
            print("initial motifs chosen")
        return seq_box

    # Use rng to make random samples/selections/numbers
    # Example: randint = rng.integer(1, 10)
    random.seed(seed)
    rng = np.random.default_rng(seed)

    seqs = utils.io_monster(mode, n_rows)

    seq_box = initialize_motifs(seqs, k)
    
    print("building frequency matrix")
    pfm = seq_box.init_pfm()                #the PFM is built according to Marcus's specs

    print("getting background frequencies")
    bg = seq_box.init_bg()
        
    pprint(bg)

    import time

    s = time.time()
    pfm = seq_box.init_pfm()
    pwm = motif_ops.build_pwm(pfm)
    e = time.time()
    print(f"{e - s} time to build pwm with seq_box")
    pprint(pwm)

    sl_seqs = seq_box.get_str_list_format_motifs()

    s = time.time()
    pfm = motif_ops.build_pfm(sl_seqs, 5)
    pwm = motif_ops.build_pwm(pfm)
    e = time.time()
    print(f"{e - s} time to build with string-list format")



    '''
    Example iteration
    '''

    for i in range(len(seq_box)):
        print(seq_box[i])
        x = utils.decode_sequence(seq_box[i])
        print(x)
        print(seqs[i])
        if x != seqs[i]:
            print("ERROR")
            input()
        print()
        input()


    '''
    -----------------------------------------------ITERATION LOOP---------------------------------------
    '''

    # converged = False
    # i = 0
    # while converged == False and i < max_iter:
    #     pfm = build_pfm()
    #     pwm = build_pwm()
    #     i += 1


GibbsMotifFinder(k=6)




# # Run the gibbs sampler:
# promoter_pfm = GibbsMotifFinder(seqs,10 )

# # Plot the final pfm that is generated: 
# seqlogo.seqlogo(seqlogo.CompletePm(pfm = promoter_pfm.T))