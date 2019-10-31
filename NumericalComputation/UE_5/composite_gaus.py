import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from numpy.polynomial.legendre import leggauss

def composite_gauss(n, L, q, f):

    # gauss points
    x_i, alpha_i = leggauss(n)

    # create intervals
    intervals = []
    intervals.append((0, q**(L-1)))
    for i in range(1,L):
        intervals.append((q**(L-i), q**(L-i-1)))

    print("intervals: ", intervals)

    f_i = []

    # integration
    for interval in intervals:
        a = interval[0]
        b = interval[1]

        # transform [a,b] to [-1,1]
        x_i_ = (b+a)/2 + x_i*(b-a)/2

        print("transformed points: ",x_i_)
        # evaluate function
        f_i.append((b-a)/2 * f(x_i_) * alpha_i)

    y = sum(sum(f_i))
    
    return y

def f(x):
    m = 4
    return np.power(x, m)

def g(x):
    return np.power(x, 0.1) * np.log(x)

y = composite_gauss(3, 5, 1/2, f)
real_y = integrate.quad(f, 0, 1)[0]
print(y - real_y)


q_list = [0.5, 0.15, 0.05]
error_list = []
true_value = -1/1.1**2
print(true_value)

L = 20
n_list = np.arange(1,21)

for q in q_list:

    temp_list = []
    for n in n_list:
        y = composite_gauss(n, L, q, g)    
        temp_list.append(abs(y-true_value))

    error_list.append(temp_list)

plt.semilogy(n_list, error_list[0], label="q = 0.5")
plt.semilogy(n_list, error_list[1], label="q = 0.15")
plt.semilogy(n_list, error_list[2], label="q = 0.05")
plt.legend()
plt.show()
