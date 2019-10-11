import numpy as np
from matplotlib import pyplot as plt

def S(N):
    numbers = np.arange(1,N+1)
    return np.sum(1.0/numbers)

b = [S(10)-np.log(10), S(100)-np.log(100), S(1000)-np.log(1000)]

A = [[1, 1.0/10, 1.0/100], [1, 1.0/100, 1.0/10000], [1, 1.0/1000, 1.0/1000000]]

a = np.linalg.solve(A, b)
print(a)
