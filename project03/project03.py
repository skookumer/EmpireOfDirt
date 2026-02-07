

import random
import numpy as np
import bamnostic as bs
import seqlogo
from pprint import pprint

#import function for building sequence motif & idenfitying seqs matching to motif
from data_readers import *
import seq_ops
import motif_ops

from joblib import Parallel, delayed

from utils import utils
from sequence_database import sequence_box
from collections import Counter, defaultdict

from scipy.special import softmax 

bases = ["A", "T", "C", "G", "N"]
encode_map = utils.init_base_encoding_map()
import time



def GibbsMotifFinder(seqs=None, k=6, n_rows=3083497, mode="norm", speed="pythonic", p_method="softmax", rtol=1e-5, atol=1e-10, max_iter=1024, seed=None, toprint=True):
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

    seqs, indptr = utils.io_monster(mode)
    motifs, midx = utils.fast_init(seqs, n_rows, indptr, k)
    seq_box = sequence_box(indptr, seqs, midx, motifs, k)

    # bg = seq_box.get_bg()


    '''
    -----------------------------------------------ITERATION LOOP---------------------------------------
    '''

    if speed == "pythonic:":

        #SLICING SUBSET
        mots = seq_box.motifs[:100]
        sample_list = []
        for j in range(len(mots)):
            sample_list.append(utils.decode_sequence(mots[j]))

        #ALTERNATE METHOD
        mots = seq_box.get_str_list_format_motifs()

        converged = False
        i = 0
        while converged == False and i < max_iter:

            mot_pick = np.random.randint(0, len(sample_list))
            sample_list.pop(mot_pick)
            pfm = motif_ops.build_pfm(sample_list, k)
            pwm = motif_ops.build_pwm(pfm)

            seq = utils.decode_sequence(seq_box[mot_pick])
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

            seq_box.midx[mot_pick] = remainder
            sample_list.insert(mot_pick, new_motif)

            if i > 0:
                if np.allclose(pwm, pwm_old, ):
                    converged = True
                    pprint(pwm-pwm_old)

            pwm_old = pwm.copy()
            i += 1
            
    elif speed == "fast":
        if toprint:
            pfm = seq_box.get_pfm()
            pprint(pfm)
        if toprint:
            print(f"beginning fast iteration")
        converged = False
        i = 0

        while converged == False and i < max_iter:
            print(f"\rIteration {i}", end="", flush=True)
            fwd_seq = seq_box.select_random_motif()
            pfm = seq_box.get_pfm(to_mask=True)
            pwm = motif_ops.build_pwm(pfm)
            rev_seq = utils.fast_complement(fwd_seq)
            fwd, rev = utils.fast_subdivide(fwd_seq, rev_seq, k)
            best_motif = utils.choose_best(fwd, rev, pwm, p_method)
            seq_box.update_motifs(best_motif)

            if i > 0:
                if np.allclose(pwm, pwm_old, rtol, atol):
                    converged = True
            pwm_old = pwm.copy()
            
            i += 1

    output = seq_box.get_str_list_format_motifs()

    if toprint:
        pfm = seq_box.get_pfm()
        print()
        pprint(pfm)
        if i == max_iter:
                print("\nalgorithm did not converge :(((")

        print(f'\nAfter {i} iterations, final motif list: {output[:20]}')
    
    return output

    
        


    

    #pprint(pfm)
    #pprint(pwm)
    #pprint(seq)
    #pprint(rev_seq)
    #pprint(scores)
    

    #pprint(sample_list)


def consensus(results, top_k=20):
    votes = Counter()
    for result in results:
        votes.update(Counter(result))
    return votes.most_common(top_k)
    

threads = 16

results = Parallel(n_jobs=threads)([delayed(GibbsMotifFinder)(k=6, speed="fast", rtol=1e-5, max_iter=6000) for _ in range(threads)])

top_k = consensus(results)
print(top_k)




# GibbsMotifFinder(k=6, speed="fast", rtol=1e-5, max_iter=6000)




# # Run the gibbs sampler:
# promoter_pfm = GibbsMotifFinder(seqs,10 )

# # Plot the final pfm that is generated: 
# seqlogo.seqlogo(seqlogo.CompletePm(pfm = promoter_pfm.T))


