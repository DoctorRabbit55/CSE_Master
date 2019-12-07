import numpy as np
import matplotlib.pyplot as plt


def f(x):
    y = np.zeros([3,1])

    y[0] = 3*x[0] - np.cos(x[1]*x[2]) - 3/2
    y[1] = 4*x[0]**2 - 625*x[1]**2 + 2*x[2] - 1
    y[2] = np.exp(-x[0]*x[1]) + 20*x[2] + 9
    return y

def df(x):
    
    y = np.matrix([[3, np.asscalar(np.sin(x[1]*x[2])*x[2]), np.asscalar(np.sin(x[1]*x[2])*x[1])],
                   [np.asscalar(8*x[0]), np.asscalar(-2*625*x[1]), 2],
                   [np.asscalar(np.exp(-x[0]*x[1])*(-x[1])), np.asscalar(np.exp(-x[0]*x[1])*(-x[0])), 20]])
    return y

N_list = np.arange(12)

error_list = []

for n in N_list:

    x0 = np.array([[1],[1],[1]])
    for i in range(n):

        f_curr = f(x0)
        df_curr = df(x0)
        delta = np.linalg.solve(df_curr, f_curr)
        x0 = x0 - delta

    error_list.append(np.linalg.norm(f(x0)))
    

plt.semilogy(N_list, error_list, label='f_error')
plt.grid()
plt.legend()
plt.show()
