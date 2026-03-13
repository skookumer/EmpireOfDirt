import math

# pairwise distance scoring functions

def hamming(seq1, seq2) -> int:
    if len(seq1) != len(seq2):
        raise ValueError("Hamming distance requires equal-length sequences")
    return sum(a != b for a, b in zip(seq1, seq2))

def p_distance(seq1: str, seq2: str) -> float:
    return hamming(seq1, seq2) / len(seq1)

def jukes_cantor(seq1: str, seq2: str, sequence_type: str = "AA") -> float:
    p = p_distance(seq1, seq2)

    if sequence_type not in ["DNA","AA"]:
        raise ValueError("sequence type should be defined (either DNA or AA) for Jukes-Cantor")

    if sequence_type == "DNA":
        if p >= 0.75:
            raise ValueError("JC69 distance undefined for p >= 0.75 for nucleotide data")
        return -0.75 * math.log(1 - (4/3) * p)
    
    # using JC for AA is not standard but here!
    return -(19/20) * math.log(1 - (4/3) * p)

def poisson(seq1: str, seq2: str) -> float:
    p = p_distance(seq1, seq2)

    if p >= 1.0:
        raise ValueError("poisson distance undefined for p >= 1 (likely it was 1, or something has gone terribly wrong).")

    return -(math.log(1.0 - p))