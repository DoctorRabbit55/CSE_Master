import numpy as np
from matplotlib import pyplot as plt

def S(N):
    numbers = np.arange(1,N+1)
    return np.sum(1.0/numbers)

def S_approx(N, a):
    return np.log(N) + a[0] + a[1]/N + a[2]/pow(N,2)

b = [S(10)-np.log(10), S(100)-np.log(100), S(1000)-np.log(1000)]

A = [[1, 1.0/10, 1.0/100], [1, 1.0/100, 1.0/10000], [1, 1.0/1000, 1.0/1000000]]

a = np.linalg.solve(A, b)

print(S_approx(10^6,a)/S(10^6))
print(S_approx(10^8,a)/S(10^8))

input('Press enter to continue..')

def S(N):
    numbers = np.power(np.arange(1,N+1),2)
    return np.sum(1.0/numbers)

b = [S(100), S(1000), S(10000)]

A = A = [[1, 1.0/100, 1.0/10000], [1, 1.0/1000, 1.0/1000000], [1, 1.0/10000, 1.0/100000000]]

a = np.linalg.solve(A, b)
print(a)
print(a[0] - pow(np.pi, 2)/6)
