from typing import Dict, List, Tuple

import numpy as np
import textdistance as td 
from pairwise_scorers import hamming, p_distance, jukes_cantor, poisson

# dispatcher that will be called in build_distance_matrix() 
# it computes pairwise distance scores using appropriate scoring model
def score_distance(seq1: str, seq2: str, scoring_model: str = "hamming", sequence_type: str = "protein")\
    -> float:
    scoring_model = scoring_model.lower()
    if scoring_model == "hamming":
        dist = hamming(seq1,seq2)
    elif scoring_model in ["p-distance", "pdistance", "p distance"]:
        dist = p_distance(seq1,seq2)
    elif scoring_model == "jukes-cantor":
        dist = jukes_cantor(seq1,seq2)
    elif scoring_model == "poisson":
        dist = poisson(seq1,seq2)
    elif scoring_model == "levenshtein":
        dist = td.levenshtein(seq1,seq2)
    elif scoring_model == "damerau-levenshtein":
        dist = td.damerau_levenshtein(seq1,seq2)
    elif scoring_model == "something else":
        pass # etc etc
    return dist

# builds D matrix by pulling helpers to compute pairwise scores 
def build_distance_matrix(sequences: Dict[str, str], scoring_model: str = "hamming", sequence_type: str = "AA")\
     -> Tuple[np.ndarray, List[str]]:
    
    ids = list(sequences.keys())
    n = len(ids)
    D = np.zeros((n, n), dtype=float) # make square n x n empty D matrix

    for i, id1 in enumerate(ids):
        for j in range(i + 1, n):
            id2 = ids[j]
            dist = score_distance(sequences[id1], sequences[id2], scoring_model)
            D[i, j] = D[j, i] = dist # mirrors value
    return D, ids

# reads fasta and gives build_distance_matrix what it wants
def fasta_to_Dict(path: str) -> Dict[str, str]:

    seqs: Dict[str, str] = {}

    with open(path, "rt", encoding="utf-8") as fh:
        current_id = None
        current_seq = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    # that means we're at the next entry
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]  # reset current id to first "word" in header
                current_seq = [] # reset current seq to empty
            else:
                current_seq.append(line)
        if current_id is not None:
            # save the last fasta entry
            seqs[current_id] = "".join(current_seq)

    return seqs

#lafayette = fasta_to_Dict("../data/lafayette_SARS_RT.fasta") 

#print(len(lafayette))
#print(lafayette.keys())
#print([len(seq) for seq in lafayette.values()]) 

# they're all 229 AA long

#D, ids = build_distance_matrix(lafayette, scoring_model="hamming")
#print(ids)
#print(D)
#row_sums_D = D.sum(axis=1)
#print(row_sums_D)

def D_to_Q(D:np.ndarray) -> np.ndarray:
    
    n = D.shape[0]
    row_sums = D.sum(axis=1) # divergence


    Q = (n - 2) * D - row_sums[:, None] - row_sums[None, :]
    # row_sums[:,None] makes a column vector of row_sums, row_sums[None, :] makes a row vector
    # numpy broadcasts column vectors across columns and row vectors across rows
    # so at Q[i,j], we've subtracted D(i,x) and D(j,x) where x is all active nodes
    # this is equivalent to "for all i, j: x_ij is x not in (i, j),  Q=sum(...)"

    np.fill_diagonal(Q, np.inf) 
    # avoid min(Q(i,j)) being the 0-distance from self to self
    # unlikely to happen but you never know

    return Q

#full_Q = D_to_Q(D)
# print(Q)

# also, since technically the only thing Q is needed for is min(Q(i,j)), we could only return Q(i,j)
# or!! even avoid ever having the full Q matrix exist entirely:

def choose_nj_pair(D: np.ndarray) -> tuple[int, int]:

    n = D.shape[0]
    row_sums = D.sum(axis=1) # divergence

    best_i, best_j = 0, 1
    best_q = np.inf

    for i in range(n):
        for j in range(i + 1, n): # this goes through the upper right triangle
            q = (n - 2) * D[i, j] - row_sums[i] - row_sums[j]
            if q < best_q:
                best_q = q
                best_i, best_j = i, j

    return best_i, best_j
    # as is, this will keep the last assessed i and j in case of ties

# best_i, best_j = choose_nj_pair(D)
# print(best_i, best_j)
# print(D[best_i,best_j])
# print(full_Q[best_i, best_j])
# print(ids[best_i])
# print(ids[best_j])

