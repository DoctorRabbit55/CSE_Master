import numpy as np
import matplotlib.pyplot as plt

def newton(x, f, df, s):
    return x - f(x, s)/df(x)
    
def f(x):
  return np.arctan(x)
  
def df(x):
  return 1/(1.0+x**2)
  
def H(x, s):
  return f(x) - (1-s) * np.arctan(4)
  
N = 3
s_list = np.arange(0,11) / 10.0
x_list = []

x_list.append(4)

for s in s_list:

  x = x_list[-1]
  for i in range(N):
    x = newton(x, H, df, s)
    
  x_list.append(x)

plt.plot(s_list, x_list[1:])
plt.show()
