import numpy as np
import pandas as pd

#for making row and col labels
def convert(string):
    list1=[]
    list1[:]=string
    list1.insert(0,0)
    return list1


def refrenceMatrix(seq1 , seq2 ,gap , mismatch1 , mismatch2, match):
    n = len(seq1)
    m = len(seq2)
    matrix = np.zeros((n,m),dtype=int)
    for i in range(n):
        for j in range(m):
            if (seq1[i] == seq2[j]):
                matrix[i][j] = match
            elif ((seq1[i] == 'A' and seq2[j] == 'G') or (seq1[i] == 'G' and seq2[j] == 'A') or
                  (seq1[i] == 'C' and seq2[j] == 'T') or (seq1[i] == 'T' and seq2[j] == 'C')):
                matrix[i][j] = mismatch1
            else:
                matrix[i][j] = mismatch2
    return matrix


def scoringMatix(seq1 , seq2 ,gap , mismatch1, mismatch2, match):
    arr = refrenceMatrix(seq1, seq2, gap, mismatch1, mismatch2, match)
    n = len(seq1)+1
    m = len(seq2)+1
    matrix = np.zeros((n,m),dtype=int)
    for i in range(n):
        matrix[i][0] = i * gap
    for j in range(m):
        matrix[0][j] = j * gap

    for i in range(1, n):
        for j in range(1, m):
            matrix[i][j] = min(matrix[i-1][j-1]+arr[i-1][j-1],
                               matrix[i-1][j]+gap,
                               matrix[i][j-1]+gap)
    return matrix

def printMatric(seq1 , seq2 ,gap , mismatch1 , mismatch2, match):
    rowlabels = convert(seq1)
    columnlabels = convert(seq2)
    mat = scoringMatix(seq1, seq2, gap, mismatch1, mismatch2, match)
    df = pd.DataFrame(mat, columns=columnlabels, index=rowlabels)
    return df


def score(seq1 , seq2 ,gap , mismatch1, mismatch2, match):
    ref = refrenceMatrix(seq1, seq2, gap, mismatch1, mismatch2, match)
    mat = scoringMatix(seq1, seq2, gap, mismatch1, mismatch2, match)

    return mat[-1][-1]


def Traceback(seq1,seq2 ,gap , mismatch1, mismatch2, match):
    ref = refrenceMatrix(seq1 , seq2 ,gap , mismatch1, mismatch2, match)
    mat = scoringMatix(seq1 , seq2 ,gap , mismatch1, mismatch2, match)

    aligned_1 = ""
    aligned_2 = ""

    ti = len(seq1)
    tj = len(seq2)

    while (ti > 0 and tj > 0):

        if (ti > 0 and tj > 0 and mat[ti][tj] == mat[ti - 1][tj - 1] + ref[ti - 1][
            tj - 1]):

            aligned_1 = seq1[ti - 1] + aligned_1
            aligned_2 = seq2[tj - 1] + aligned_2

            ti = ti - 1
            tj = tj - 1

        elif (ti > 0 and mat[ti][tj] == mat[ti - 1][tj] + gap):
            aligned_1 = seq1[ti - 1] + aligned_1
            aligned_2 = "-" + aligned_2

            ti = ti - 1
        else:
            aligned_1 = "-" + aligned_1
            aligned_2 = seq2[tj - 1] + aligned_2

            tj = tj - 1
    return aligned_1 ,aligned_2





def main():
    #dynamic all that could be taken from user
    seq1 = 'TACGTCAGC'
    seq2 = 'TATGTCATGC'
    gap = 8
    mismatch1 = 2
    mismatch2 = 4
    match = 0

    newLine = '\n'
    result = score(seq1 , seq2 ,gap , mismatch1 , mismatch2, match)
    mat = printMatric(seq1 , seq2 ,gap , mismatch1 , mismatch2, match)
    out = Traceback(seq1,seq2 ,gap , mismatch1, mismatch2, match)
    print(f"Optimal global alignment value is {result} ,{newLine}and the matrix is {newLine}{mat}{newLine}{out[0]}{newLine}{out[1]}")


import unittest


class Test(unittest.TestCase):
    def test(self):
        seq1 = 'TACGTCAGC'
        seq2 = 'TATGTCATGC'
        gap = 8
        mismatch1 = 2
        mismatch2 = 4
        match = 0
        result = score(seq1, seq2, gap, mismatch1, mismatch2, match)
        self.assertEqual(result,10)


if __name__ == "__main__":
    main()
    unittest.main()

