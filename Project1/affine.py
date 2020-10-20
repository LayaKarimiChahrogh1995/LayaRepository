import numpy as np
from pam_n import pam_n


subs_matrix = pam_n(250)

def seq_mapping():
    index_mapping = {}
    index_mapping ['A'] = 0
    index_mapping ['R'] = 1
    index_mapping ['N'] = 2
    index_mapping ['D'] = 3
    index_mapping ['C'] = 4
    index_mapping ['Q'] = 5
    index_mapping ['E'] = 6
    index_mapping ['G'] = 7
    index_mapping ['H'] = 8
    index_mapping ['I'] = 9
    index_mapping ['L'] = 10
    index_mapping ['K'] = 11
    index_mapping ['M'] = 12
    index_mapping ['F'] = 13
    index_mapping ['P'] = 14
    index_mapping ['S'] = 15
    index_mapping ['T'] = 16
    index_mapping ['W'] = 17
    index_mapping ['Y'] = 18
    index_mapping ['V'] = 19
    return index_mapping

def match_score(alpha, beta, index_mapping):
    return subs_matrix[index_mapping[alpha], index_mapping[beta]]

seq1 = "ATTACA"
seq2 = "ATGCT"


def affine_align(seq1, seq2):
    gap_penalty = int(input('Enter your gap penalty :\n'))
    gap_start_penalty = int(input('Enter your gap start penalty score :\n'))
    gap_extend_penalty = int(input('Enter your gap extend penalty score :\n'))
    m = len(seq1)
    n = len(seq2)
    M = np.zeros((len(seq1) + 1, len(seq2) + 1))
    X = np.zeros((len(seq1) + 1, len(seq2) + 1))
    Y = np.zeros((len(seq1) + 1, len(seq2) + 1))

    for i in range(1, m + 1):
        M[i][0] = float('-inf')
        X[i][0] = gap_start_penalty + i * gap_penalty
        Y[i][0] = float('-inf')

    for i in range(1, n + 1):
        M[0][i] = float('-inf')
        X[0][i] = gap_start_penalty + i * gap_penalty
        Y[0][i] = float('-inf')
    index_mapping = seq_mapping()
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            M[i][j] = match_score(seq1[i - 1], seq2[j - 1], index_mapping) + max(M[i - 1][j - 1], X[i - 1][j - 1],
                                                                                 Y[i - 1][j - 1])
            X[i][j] = max(gap_start_penalty + gap_extend_penalty + M[i][j - 1],
                          gap_extend_penalty + X[i][j - 1],
                          gap_start_penalty + gap_extend_penalty + Y[i][j - 1])
            Y[i][j] = max(gap_start_penalty + gap_extend_penalty + M[i - 1][j],
                          gap_start_penalty + gap_extend_penalty + X[i - 1][j],
                          gap_extend_penalty + Y[i - 1][j])
    opt = max(M[m][n], X[m][n], Y[m][n])
    print('Optimal score', opt)
    print('M matrix \n:', M)
    print('Y matrix: \n', Y)
    print('X matrix:\n', X)

    index_mapping = seq_mapping()
    align1 = ''
    align2 = ''
    i = len(seq1)
    j = len(seq2)

    while (i > 0 or j > 0):
        score_current = opt
        score_diagonal = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1])
        score_up = max(M[i][j - 1], X[i][j - 1], Y[i][j - 1])
        score_left = max(M[i - 1][j], X[i - 1][j], Y[i - 1][j])

        if score_current == score_diagonal + match_score(seq2[j - 1], seq1[i - 1], index_mapping):
            align1 += seq2[j - 1]
            align2 += seq1[i - 1]
            i -= 1
            j -= 1
        elif score_current == score_up + gap_penalty or score_up + gap_start_penalty or score_up + gap_extend_penalty:
            align1 += seq2[j - 1]
            align2 += '-'
            j -= 1
        elif score_current == score_left + gap_penalty or score_up + gap_start_penalty or score_up + gap_extend_penalty:
            align1 += '-'
            align2 += seq1[i - 1]
            i -= 1
        # Finish tracing up to the top left cell
        while j > 0:
            align1 += seq2[j - 1]
            align2 += '-'
            j -= 1
        while i > 0:
            align1 += '-'
            align2 += seq1[i - 1]
            i -= 1
        align1 = align1[::-1]
        align2 = align2[::-1]
        return (align1, align2)



print(affine_align(seq1,seq2))
