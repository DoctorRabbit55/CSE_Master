import numpy as np
import matplotlib.pyplot as plt

def newton(x, f, df):
    return x - f(x)/df(x)


def f1(x):
    return np.power(x, 2)

def df1(x):
    return 2*x

def f2(x):
    return np.exp(x)-2

def df2(x):
    return np.exp(x)

def f3(x):
    return np.power(np.abs(x), 3.0/2.0)

def df3(x):
    if x == 0:
        return 0

    return 3*x/(2*np.power(np.abs(x), 1.0/2.0))

def f4(x):
    return 1/x - 1

def df4(x):
    return -1/np.power(x,2)


N_list = np.arange(1, 15, 2).astype(np.int16)
error_list_1 = []

for N in N_list:

    x0 = 0.5
    for i in range(N):
        x0 = newton(x0, f1, df1)

    error_list_1.append(abs(x0))

error_list_2 = []
for N in N_list:

    x0 = 0.5
    for i in range(N):
        x0 = newton(x0, f2, df2)

    error_list_2.append(abs(np.log(2)-x0))

error_list_3 = []
for N in N_list:

    x0 = 0.5
    for i in range(N):
        x0 = newton(x0, f3, df3)

    error_list_3.append(abs(x0))

error_list_4 = []
for N in N_list:
    x0 = 2.1
    for i in range(N):
        x0 = newton(x0, f4, df4)

    error_list_4.append(abs(1-x0))

plt.subplot(221)
plt.semilogy(N_list, error_list_1, label='f1')
plt.semilogy(N_list, 1.0/np.power(2, N_list), label='O(N^2)')
plt.grid()
plt.legend()
plt.subplot(222)
plt.semilogy(N_list, error_list_2, label='f2')
plt.semilogy(N_list, 1.0/np.power(2, N_list), label='O(N^2)')
plt.grid()
plt.legend()
plt.subplot(223)
plt.semilogy(N_list, error_list_3, label='f3')
plt.semilogy(N_list, 1.0/np.power(2, N_list), label='O(N^2)')
plt.grid()
plt.legend()
plt.subplot(224)
plt.semilogy(N_list, error_list_4, label='f4')
plt.semilogy(N_list, 1.0/np.power(2, N_list), label='O(N^2)')
plt.grid()
plt.legend()
plt.show()
    
