

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
import time



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

    # Use rng to make random samples/selections/numbers
    # Example: randint = rng.integer(1, 10)
    random.seed(seed)
    rng = np.random.default_rng(seed)

    seqs, indptr = utils.io_monster(mode, n_rows)
    motifs, midx = utils.fast_init(seqs, n_rows, indptr, k)
    seq_box = sequence_box(indptr, seqs, midx, motifs, k)

    bg = seq_box.init_bg()
    pfm = seq_box.init_pfm()
    pwm = motif_ops.build_pwm(pfm)

    pprint(pfm)


    '''
    Example iteration
    '''

    for i in range(len(seq_box)):
        print(seq_box[i])
        x = utils.decode_sequence(seq_box[i])
        print(x)
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