import numpy as np
import pandas



var1 = np.array([[1,2,3],[2,3,4]])
var2 = np.zeros((100,100))

seq1 = 'ATCT'
seq2 = 'TCTA'


seqMat = np.zeros((len(seq1), len(seq2)), dtype=np.int)


def func(i):
    return i%4+1

npfunc = np.fromfunction(func, (10,), dtype=np.int)
I = numpy.where(np='xxx')
