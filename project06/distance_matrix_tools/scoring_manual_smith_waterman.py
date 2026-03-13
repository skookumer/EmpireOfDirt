import numpy as np
import textdistance

# manual smith-waterman
# this returns similarity/alignment scores which can't be used as is in a distance matrix
# but it's here if we want it 

# helper: the cell-wise smith-waterman operation, returns (score, direction)
def _cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap):

    # diagonal travel (match/mismatch)
    if seq1[i - 1] == seq2[j - 1]:
        diag_score = matrix[i - 1, j - 1] + match # match
    else:
        diag_score = matrix[i - 1, j - 1] + mismatch # mismatch

    # seq1[i-1] is a aligned with gap in seq2
    up_score = matrix[i - 1, j] + gap

    # seq2[j-1] is a aligned with gap in seq1
    left_score = matrix[i, j - 1] + gap

    # score
    score = max(
        0, # all options are bad, reset to 0
        diag_score, 
        up_score, 
        left_score
        )

    # store which move gave the best score
    if score == 0:
        direction = "end"
    elif score == diag_score:
        direction = "diag"
    elif score == up_score:
        direction = "up"
    else:
        direction = "left"

    return score, direction

# helper: get alignment by backtracking from a cell in built smith-waterman matrix
def _traceback(seq1, seq2, traceback_matrix, maximum_position):
    aligned_seq1 = []
    aligned_seq2 = []

    i, j = maximum_position

    while traceback_matrix[i, j] != "end":
        direction = traceback_matrix[i, j]

        if direction == "diag":
            # best score came from diag
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == "up":
            # best score came from up = gap in seq2
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        elif direction == "left":
            # best score came from left = gap in seq1
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    aligned_seq1.reverse()
    aligned_seq2.reverse()

    # join list into a string
    return "".join(aligned_seq1), "".join(aligned_seq2)

# helper: get smith-waterman score for two sequences, defaults to only returning a score
def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-1, score_only=True):
    
    # prep "empty" score matrix and traceback matrix to be updated
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score_matrix = np.zeros((rows, cols), dtype=int) # all 0 scores
    traceback_matrix = np.full((rows, cols), "end", dtype=object) # all "end" directions

    # prep max score and max position to be updated
    max_score = 0
    max_position = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):

            # cal_score computes score and direction
            score, direction = _cal_score(score_matrix, seq1, seq2, i, j, match, mismatch, gap)

            # store outputs
            score_matrix[i, j] = score
            traceback_matrix[i, j] = direction

            # update max_score and max_position if new max score found
            if score > max_score:
                max_score = score
                max_position = (i, j)

    if score_only:
        return max_score

    # get the aligned sequences, if you want
    aligned_seq1, aligned_seq2 = _traceback(seq1, seq2, traceback_matrix, max_position)
    return (max_score, aligned_seq1, aligned_seq2, score_matrix)