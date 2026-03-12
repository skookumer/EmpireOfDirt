# Import package
import numpy as np  # array for dp

# Direction codes
END, DIAG, UP, LEFT = 0, 1, 2, 3


def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap):
    """
    Compute score for cell (i,j), return score and traceback direction
    Parameters:
        matrix (np.array): numpy 2d scoring matrix
        seq1 (str): vertical sequence (rows)
        seq2 (str): horizontal sequence (columns)
        i, j (int): current cell indices
    Returns:
        score (float): max score in matrix
        traceback (int): traceback direction
    """

    # diag-> match/mismatch
    # Add match reward if same, mismatch if not but diag still max
    diag = matrix[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)

    # up-> gap in seq2
    up = matrix[i-1, j] + gap

    # left-> gap in seq1
    left = matrix[i, j-1] + gap

    # local alignment takes max score of zero or scores from 3 directions
    score = max(0, diag, up, left)

    # determine traceback direction
    # *** FIXED: DIAG no longer wins ties when mismatch = gap ***
    if score == 0:
        direction = END  # Always stop alignment once zero reached

    # If UP and LEFT tie and both equal score, choose stronger neighbor
    elif score == up and score == left:
        if matrix[i-1, j] >= matrix[i, j-1]:
            direction = UP
        else:
            direction = LEFT

    # If UP is the unique max
    elif score == up:
        direction = UP

    # If LEFT is the unique max
    elif score == left:
        direction = LEFT

    # DIAG chosen ONLY if strictly the max (match case)
    else:
        direction = DIAG  # can be match or mismatch

    return score, direction


def traceback(seq1, seq2, score_matrix, traceback_matrix, start_pos):
    """
    Reconstruct optimal local alignment using Smith–Waterman
    traceback rule: traceback begins at highest-scoring cell,
    proceeds until cell containing score 0 reached
    Parameters:
        seq1 (str): Vertical sequence (rows)
        seq2 (str): Horizontal sequence (columns)
        score_matrix (np.array): Scoring matrix produced during alignment
        traceback_matrix (np.array): Matrix of direction codes (END, DIAG, UP, LEFT) for traceback
        start_pos (tuple[int, int]): Highest scoring cell coordinates in scoring matrix
    Returns:
        aligned1 (str): Aligned substring of seq1, including gap characters
        aligned2 (str): Aligned substring of seq2, including gap characters
    """

    # Initialize empty lists (not strings) for aligned seqs
    aligned1 = []
    aligned2 = []

    # Set starting coordinates for traceback as max score in matrix
    i, j = start_pos

    # Traceback continues until zero score encountered
    while score_matrix[i, j] != 0:

        # direction code for this cell in traceback matrix
        move = traceback_matrix[i, j]

        # Diagonal move-> match or mismatch
        if move == DIAG:
            aligned1.append(seq1[i - 1])  # Consume character from seq1 for alignment
            aligned2.append(seq2[j - 1])  # Consume character from seq2 for alignment
            i -= 1
            j -= 1

        # Up move
        elif move == UP:
            aligned1.append(seq1[i - 1])  # seq1 gets real character
            aligned2.append('-')  # Insert gap in seq2
            i -= 1  # move traceback pointer up

        # Left move
        elif move == LEFT:
            aligned1.append('-')  # Insert gap in seq1
            aligned2.append(seq2[j - 1])  # seq2 gets real character
            j -= 1  # move traceback pointer left

        # Fallback termination (rare)-> stop if direction is END
        else:
            break

    # Reverse lists-> traceback builds alignment backwards
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))

    return aligned1, aligned2


def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1):
    """
    Smith–Waterman local alignment using NumPy array as matrix
    Parameters:
        seq1: (str): string sequence of bases to align
        seq2: (str): string sequence of bases to align
        match: (int): score reward for matching alignments
        mismatch: (int): score penalty for mismatching alignments
        gap (int): score penalty for sequence of unequal lengths
    Returns:
        aligned seq1: (str): string sequence of bases from seq1 with
                             highest alignment score
        aligned seq2: (str): string sequence of bases from seq2 with
                             highest alignment score
        score matrix: (array): complete alignment score matrix
    """

    # matrix dimensions: (len(seq1)+1) x (len(seq2)+1) -> add 1 for initialization row/col
    rows, cols = len(seq1) + 1, len(seq2) + 1

    # Initialize scoring and traceback matrices with zeros
    score_matrix = np.zeros((rows, cols), dtype=int)
    traceback_matrix = np.zeros((rows, cols), dtype=int)

    # Initialize max score and matrix position variables
    max_score = 0  # maximum score in matrix
    max_pos = (0, 0)  # coordinates of max position in matrix

    # Construct DP matrix
    for i in range(1, rows):
        for j in range(1, cols):
            score, direction = cal_score(score_matrix, seq1, seq2, i, j,
                                         match, mismatch, gap)
            # Calculate score
            score_matrix[i, j] = score
            # Determine direction
            traceback_matrix[i, j] = direction

            # Track global maximum for local alignment
            if score > max_score:
                max_score = score
                max_pos = (i, j)

    # Traceback from highest scoring cell
    aligned1, aligned2 = traceback(seq1, seq2, score_matrix, traceback_matrix, max_pos)

    return aligned1, aligned2, score_matrix


if __name__ == "__main__":

    # Perfect match alignment test
    # Sanity check-> ensures full diagonal match & no zero resets
    # seq1 = "ATGCGTACGTTAGCTAGCTAGGCTAACGTTAGC"
    # seq2 = "ATGCGTACGTTAGCTAGCTAGGCTAACGTTAGC"
    # print("Expected Result: Full Alignment")

    # Long sequences with one local island
    # tests ability to isolate a long internal match
    # ensures flanking noise doesn’t inflate alignment
    # ensures traceback starts/ends correctly
    # seq1 = "GGGGGGGGGGACGTACGTACGTACGTGGGGGGGGGG"
    # seq2 = "TTTTTTTTTTACGTACGTACGTACGTTTTTTTTTTT"
    # print("Expected Result: Local Alignment ACGTACGTACGTACGT")

    # Multiple competing local maxima -> 3 equally strong local alignments
    # seq1 = "ACGTACGTTTTTTTTTTTTACGTACGTTTTTTTTTTTTACGTACGT"
    # seq2 = "ACGTACGTGGGGGGGGGGGGACGTACGTGGGGGGGGGGGGACGTACGT"
    # print("Expected Result: single local alignment maximum ACGTACGT)")

    # Zero reset deep in matrix ->
    # long mismatch regions force many zero resets
    # local match in the middle
    # ensures traceback stops at correct zero boundary
    # seq1 = "TTTTTTTTTTTTTTTTACGTTGCAACGTGGGGGGGGGGGG"
    # seq2 = "GGGGGGGGGGGGGGGGACGTTGCAACGTAAAAAAAAAAAA"
    # print("Expected Result: ACGTTGCAACGT only")

    # Long seqs with only one shared base ->
    # tests minimal alignment
    # ensures no accidental extension
    # checks correct traceback termination
    # seq1 = "GGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
    # seq2 = "TTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    # print("Expected Result: A only")

    # Biological-style conserved motif (long)
    # tests mixed match/mismatch
    # tests traceback through a conserved region
    # tests that local alignment finds the motif, not the ends
    seq1 = "ATGCGTACCTGACCTGACCTGACCTGACCTGA"
    seq2 = "TTACGTACCTGACCTGACCTGACCTGACCGGG"
    print("Expected Result: ACGTACCTGACCTGACCTGACCTGACCT")

    # Internal inversion test ->
    # Program must not align reverse complements
    # Ensures local alignment doesn’t incorrectly match inverted regions
    # seq1 = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    # seq2 = "TGCATGCATGCATGCATGCATGCATGCATGCA"
    # print("Expected Result: Very short or zero alignment")

    # Long repetitive sequences/ambiguous diagonals
    # Tests DP behavior on repeats
    # Tests diagonal ambiguity
    # Tests tie‑breaking in repetitive regions
    # seq1 = "ATATATATATATATATATATATATATATATATATAT"
    # seq2 = "ATATATATATATATATATATATATATATATATATAT"
    # print("Expected Result: Full ATATAT... alignment")

    # Long gap inside high scoring region
    # Tests gap handling inside a conserved region
    # Ensures traceback doesn’t break long gaps into fragments
    # Exposes gap‑penalty mis‑implementations
    # seq1 = "ACGTAAAAAACGTAAAAAACGTAAAAAACGT"
    # seq2 = "ACGTACGTACGTACGTACGTACGTACGT"
    # print("Expected Result: ACGT----ACGT----ACGT----ACGT")

    # Multiple local alignments with different scores
    # Tests that the algorithm picks the highest scoring island
    # Not the longest, first or leftmost
    # seq1 = "AAAAACGTAAAAACGTAAAAACGTAAAAACGTAAAAA"
    # seq2 = "TTTTACGTGGGGACGTGGGGACGTGGGGACGTTTTTT"
    # print("Expected Result: Highest scoring ACGT island only")

    aligned1, aligned2, score_matrix = smith_waterman(seq1, seq2)

    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print("Aligned seq1:", aligned1)
    print("Aligned seq2:", aligned2)
    np.set_printoptions(threshold=np.inf)
    print("Score matrix:\n", score_matrix)
