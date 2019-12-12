import numpy as np
import matplotlib.pyplot as plt
import math


def f(x,h):
    
    y = np.zeros((len(x), 1), dtype=np.float64)
    
    for i in range(len(x)):
        if(i==0):
            y[i] = (-x[i+1] + 2*x[i]) / h**2 + x[i]**3 - 1
        if(i==len(x)-1):
            y[i] = (2*x[i] - x[i-1])/h**2 + x[i]**3 - 1
        else:
            y[i] = (-x[i+1] + 2*x[i] - x[i-1])/h**2 + x[i]**3 - 1
    return y

def df(x,h):

    y = np.zeros((len(x), len(x)), dtype=np.float64)
    for i in range(len(x)):
        
        if (i == 0):
            
            y[i][i] = 2/h**2+3*x[i]**2
            y[i][i+1] = -1/h**2
            
        if (i == len(x)-1):
            
            y[i][i-1] = -1/h**2
            y[i][i] = 2/h**2+3*x[i]**2

        else:
            
            y[i][i-1] = -1/h**2
            y[i][i] = 2/h**2+3*x[i]**2
            y[i][i+1] = -1/h**2
        
    return y



## Code ##
    
N = 10
h = 4.0
iterations = 20
error_list = []

x0 = np.random.rand(N, 1)

for it in range(iterations):

    print(x0)
    f_curr = f(x0, h)
    df_curr = df(x0, h)
    delta = np.linalg.solve(df_curr, f_curr)
    x0 = x0 - delta
    error_list.append(np.linalg.norm(delta))
    
plt.semilogy(range(iterations), error_list, label='error')
plt.grid()
plt.xlabel("iteration number")
plt.ylabel("error")
plt.legend()
plt.show()
