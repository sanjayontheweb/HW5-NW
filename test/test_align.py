# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

#Assisted by Chat-Gpt-4o

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10, -1)
    score, align1, align2 = nw.align(seq1, seq2)

    align_matrix = np.array([
       [  0., -np.inf, -np.inf, -np.inf],
       [-np.inf,   5., -11., -13.],
       [-np.inf, -12.,   4.,  -8.],
       [-np.inf, -12.,  -1.,   5.],
       [-np.inf, -14.,  -6.,   4.]])
    
    # Test final alignment
    assert np.allclose(nw._align_matrix, align_matrix)
    assert score == 4.0

    
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", -10, -1)
    score, align1, align2 = nw.align(seq3, seq4)

    assert score == 17
    assert align1 == 'MAVHQLIRRP'
    assert align2 == 'M---QLIRHP'



