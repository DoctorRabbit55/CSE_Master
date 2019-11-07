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

n_list = np.arange(1,21)
L_list = n_list

for q in q_list:

    temp_list = []
    for i in range(len(n_list)):
        y = composite_gauss(n_list[i], L_list[i], q, g)    
        temp_list.append(abs(y-true_value))

    error_list.append(temp_list)

plt.semilogy(n_list, error_list[0], label="q = 0.5")
plt.semilogy(n_list, error_list[1], label="q = 0.15")
plt.semilogy(n_list, error_list[2], label="q = 0.05")
plt.legend()
plt.show()

error_list = np.asarray(error_list)
error_list[1] = error_list[1] + 1e-16

error_list_log = []
error_list_log.append(np.log(error_list[0]))
error_list_log.append(np.log(error_list[1]))
error_list_log.append(np.log(error_list[2]))

fit_1 = np.polyfit(n_list, error_list_log[0], 1)
fit_2 = np.polyfit(n_list, error_list_log[1], 1)
fit_3 = np.polyfit(n_list, error_list_log[2], 1)

print("q=0.5:",fit_1)
print("q=0.15:",fit_2)
print("q=0.05:",fit_3)

fit_1 = np.poly1d(fit_1)
fit_2 = np.poly1d(fit_2)
fit_3 = np.poly1d(fit_3)

plt.semilogy(n_list, np.exp(fit_1(n_list)), color='C0')
plt.semilogy(n_list, np.exp(fit_2(n_list)), color='C1')
plt.semilogy(n_list, np.exp(fit_3(n_list)), color='C2')
plt.semilogy(n_list, error_list[0], label="q = 0.5", linestyle=":", color='C0')
plt.semilogy(n_list, error_list[1], label="q = 0.15", linestyle=":", color='C1')
plt.semilogy(n_list, error_list[2], label="q = 0.05", linestyle=":", color='C2')
plt.legend()
plt.show()
