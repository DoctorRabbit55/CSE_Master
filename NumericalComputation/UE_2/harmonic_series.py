import numpy as np
import time
from matplotlib import pyplot as plt

def S(N):
    numbers = np.arange(1,N+1)
    return np.sum(1.0/numbers)

def S_approx(N, a):
    return np.log(N) + a[0] + a[1]/N + a[2]/pow(N,2)

b = [S(10)-np.log(10), S(100)-np.log(100), S(1000)-np.log(1000)]

A = [[1, 1.0/10, 1.0/100], [1, 1.0/100, 1.0/10000], [1, 1.0/1000, 1.0/1000000]]

a = np.linalg.solve(A, b)

time1 = time.time()
S_ap = S_approx(10**6, a)
time2 = time.time()
S_real = S(10**6)
time3 = time.time()

print("time of approximation:",time2-time1)
print("time of real sum:",time3-time2)
print("error of approximation:",S_ap/S_real)

input('Press enter to continue..')

time1 = time.time()
S_ap = S_approx(10**8, a)
time2 = time.time()
S_real = S(10**8)
time3 = time.time()

print("time of approximation:",time2-time1)
print("time of real sum:",time3-time2)
print("error of approximation:",S_ap/S_real)

input('Press enter to continue..')

def S(N):
    numbers = np.power(np.arange(1,N+1),2)
    return np.sum(1.0/numbers)

b = [S(100), S(1000), S(10000)]

A = A = [[1, 1.0/100, 1.0/10000], [1, 1.0/1000, 1.0/1000000], [1, 1.0/10000, 1.0/100000000]]

a = np.linalg.solve(A, b)
print(a)
print(a[0] - pow(np.pi, 2)/6)
