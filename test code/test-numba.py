# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 20:21:52 2014

@author: pedro
"""

from __future__ import division
from numba import jit
from numpy import arange
import time

# jit decorator tells Numba to compile this function.
# The argument types will be inferred by Numba when function is called.
@jit
def sum2d(arr):
    M, N = arr.shape
    result = 0.0
    for i in range(M):
        for j in range(N):
            result += arr[i,j]
    return result

ini = time.time()
a = arange(9000000).reshape(3000,3000)
ini2 = time.time()
res = sum2d(a)
print time.time() - ini2
print res