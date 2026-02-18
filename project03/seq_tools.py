import numpy as np
import random
from scipy.special import softmax


def read_file(filename):
    '''Takes file name and returns list of sequences found in file
    Params: filename (string) - name of file to be opened
    Return: seqs (list) - list of strings of nucleotides
    '''

    seqs = []
    with open(filename) as file:
        for line in file:
            seqs.append(line.strip().strip("N"))

    return seqs


def build_motif_starts(peaks, k):
    '''Build the initial motif start sites for the initialization
    Params: peaks (list) - list of strings of nucleotides
            k (int) - k-mer length
    Return: starts (list) - list of ints representing k-mer start indices
    '''
    starts = []
    for seq in peaks:
        # makes sure to only choose a possible start site up to len-k
        starts.append(random.randrange(len(seq)-k))

    return starts


def get_kmer_list(seqs, starts, k):
    '''Take in sequences and start locations and return all the chosen kmers
    Params: seqs (list) - list of strings of nucleotides
            starts (list) - list of ints representing k-mer start indices
    Return: kmer_list (list) - list of strings representing k-mers
    '''
    kmer_list = []
    i = 0
    for seq, idx in zip(seqs, starts):
        kmer_list.append(seq[idx:idx+k])
        i += 1

    return kmer_list


def change_pfm(sequence, k, operation, pfm):
    '''Either add to or remove one k-mer worth of sequence information from pfm
    Params: sequence: (string) motif to add or remove
            k: (int) length of kmer
            pfm (np.array) - position frequency matrix
            operation: (str) "add" or "sub" from pfm
    '''
    # does not require return value becasue np.add and np.subtract directly affect the np array.
    
    base_to_index = np.zeros(256, dtype=np.int8)
    base_to_index[ord('A')] = 0
    base_to_index[ord('C')] = 1
    base_to_index[ord('G')] = 2
    base_to_index[ord('T')] = 3

    seq_array = np.frombuffer(sequence.encode(), dtype=np.int8)
    indices = base_to_index[seq_array]
    if operation == "add":
        np.add.at(pfm, (indices, np.arange(k)), 1)
    elif operation == "sub":
        np.subtract.at(pfm, (indices, np.arange(k)), 1)
