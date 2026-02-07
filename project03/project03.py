

import random
import numpy as np
import bamnostic as bs
import seqlogo
from pprint import pprint

#import function for building sequence motif & idenfitying seqs matching to motif
from data_readers import *
import seq_ops
import motif_ops

#from joblib import parallel, delayed

from utils import utils
from sequence_database import sequence_box
from collections import Counter, defaultdict

from scipy.special import softmax 

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

    #bg = seq_box.init_bg()
    #pfm = seq_box.init_pfm()
    #pwm = motif_ops.build_pwm(pfm)

    #pprint(pfm)


    '''
    Example iteration
    '''

    #for i in range(len(seq_box)):
     #  x = utils.decode_sequence(seq_box[i])
      #  print(x)
       # input()

    bg = seq_box.get_bg()
    pfm = seq_box.get_pfm()
    pwm = motif_ops.build_pwm(pfm)


    '''
    -----------------------------------------------ITERATION LOOP---------------------------------------
    '''

    seqs = seq_box.motifs[:100]
    sample_list = []
    for j in range(len(seqs)):
        sample_list.append(utils.decode_sequence(seqs[j]))

    converged = False
    i = 0
    while converged == False and i < max_iter:

        seq_pick = np.random.randint(0, len(sample_list))
        sample_list.pop(seq_pick)
        pfm = motif_ops.build_pfm(sample_list, k)
        pwm = motif_ops.build_pwm(pfm)

        seq = utils.decode_sequence(seq_box[seq_pick])
        rev_seq = seq_ops.reverse_complement(seq)

        scores = []
        for x in range(len(seq) - k):
            score = motif_ops.score_kmer(seq[x:x+k], pwm)
            scores.append(score)
        
        for x in range(len(seq) - k):
            score = motif_ops.score_kmer(rev_seq[x:x+k], pwm)
            scores.append(score)
        
        prob = softmax(scores)
        #print(prob)
        idx = np.random.choice(np.arange(len(scores)), p=prob)

        quotient, remainder = divmod(idx, len(seq) - k - 1)
        if quotient == 0:
            new_motif = seq[remainder:remainder+k]
        else:
            new_motif = rev_seq[remainder:remainder+k]

        seq_box.midx[seq_pick] = remainder
        sample_list.insert(seq_pick, new_motif)

        if i > 0:
            if np.allclose(pwm, pwm_old):
                converged = True
                pprint(pwm-pwm_old)

        pwm_old = pwm.copy()
        i += 1

    print(f'After {i} iterations, final motif list: {sample_list}')
    
        


    

    #pprint(pfm)
    #pprint(pwm)
    #pprint(seq)
    #pprint(rev_seq)
    #pprint(scores)
    

    #pprint(sample_list)

GibbsMotifFinder(k=6)




# # Run the gibbs sampler:
# promoter_pfm = GibbsMotifFinder(seqs,10 )

# # Plot the final pfm that is generated: 
# seqlogo.seqlogo(seqlogo.CompletePm(pfm = promoter_pfm.T))


