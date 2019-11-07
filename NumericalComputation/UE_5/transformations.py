import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from numpy.polynomial.legendre import leggauss

def composite_gauss(n, L, q, f, bounds):

    # gauss points
    x_i, alpha_i = leggauss(n)

    # create intervals
    intervals = []
    intervals.append((bounds[0], bounds[1]*q**(L-1)))
    for i in range(1,L):
        intervals.append((bounds[1]*q**(L-i), bounds[1]*q**(L-i-1)))

    #print("intervals: ", intervals)

    f_i = []

    # integration
    for interval in intervals:
        a = interval[0]
        b = interval[1]

        # transform [a,b] to [-1,1]
        x_i_ = (b+a)/2 + x_i*(b-a)/2

        #print("transformed points: ",x_i_)
        # evaluate function
        f_i.append((b-a)/2 * sum(f(x_i_)*alpha_i))


    y = sum(f_i)
    
    return y


def f(x):
    return np.log(1/x)*1/np.power(1/x, np.pi-2)

def g(x):
    return np.log(np.exp(x))*1/np.power(np.exp(x), np.pi)*np.exp(x)

q = 0.15
error_list_1 = []
error_list_2 = []
true_value = 1/(np.pi**2 - 2*np.pi +1)

n_list = np.arange(1,21)
L_list = n_list

for i in range(len(n_list)):
    y = composite_gauss(n_list[i], L_list[i], q, f, [0, 1])    
    error_list_1.append(abs(y-true_value))
    
    y = composite_gauss(n_list[i], L_list[i], q, g, [0, 100])    
    error_list_2.append(abs(y-true_value))


plt.semilogy(n_list, error_list_1, label="y = 1/x")
plt.semilogy(n_list, error_list_2, label="y = e^x")
plt.legend()
plt.show()


