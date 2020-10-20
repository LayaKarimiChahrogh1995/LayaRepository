import numpy as np
from pam_n import pam_n

subs_matrix = pam_n(250)

# Use these values to calculate scores
#gap_penalty = int(input('Enter your gap penalty :\n'))
#match_award = int(input('Enter your match score :\n'))
#mismatch_penalty = int(input('Enter your mismatch score :\n'))

seq1 = "IETVSY"
seq2 = "SLVDM"



def zeros(rows, cols):
    retval = []
    for x in range(rows):
        retval.append([])
        for y in range(cols):
            retval[-1].append(0)
    return retval

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

#determining the score between any two bases of alignment
def match_score(alpha, beta, index_mapping):
    return subs_matrix[index_mapping[alpha], index_mapping[beta]]

# filling out the matrix of scores
def semi_global_alignment(seq1, seq2):
    gap_penalty = int(input('Enter your gap penalty :\n'))
 # Store length of two sequences
    n = len(seq1)
    m = len(seq2)

    score = zeros(m+1, n+1)

    # Calculating score table

    for i in range(0, m + 1):
        score[i][0] = 0

    for j in range(0, n + 1):
        score[0][j] = 0

    index_mapping = seq_mapping()
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[j-1], seq2[i-1], index_mapping)
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)
    semi_global = np.array(score)
    print(semi_global)
    opt_score = []
    best = 0
    optloc = (0,0)

    for i in range (1, m + 1):
        values = score[i][m+1]
        if score[i][m+1] >=best:
            best = score[i][m+1]
            optloc = (i,j)
        opt_score.append(values)

    for j in range (1, n + 1):
        values = score[n-1][j]
        if score [n-1][j] >= best:
            best = score[n-1][j]
            optloc = (i,j)
        opt_score.append(values)
    print(optloc)
    print('the opt score value is:', max(opt_score) )
    print('the values from final row and column',opt_score)
    print(optloc)
    index_max = np.array(optloc)


    align1 = ""
    align2 = ""
    i = index_max[0]
    j = index_max[1]
    while (i > 0 and j > 0) and score[i][j] > 0: # end touching the top or the left edge
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]

        if score_current == score_diagonal + match_score(seq1[j-1],seq2[i-1], index_mapping):
            align1 += seq1[j-1]
            align2 += seq2[i-1]
            i -= 1
            j -= 1
        elif score_current == score_up + gap_penalty:
            align1 += seq1[j-1]
            align2 += '-'
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += '-'
            align2 += seq2[i-1]
            i -= 1

    while j > 0:
        align1 += seq1[j-1]
        align2 += '-'
        j -= 1
    while i > 0:
        align1 += '-'
        align2 += seq2[i-1]
        i -= 1

    align1 = align1[::-1]
    align2 = align2[::-1]
    return (align1,align2)

result1, result2 = semi_global_alignment(seq1,seq2)
print("Semi global alignment seq is: \n" + result1 + '\n' + result2)
#semi_global_alignment(seq1, seq2)
