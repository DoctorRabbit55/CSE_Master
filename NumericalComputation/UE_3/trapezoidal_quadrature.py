import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def f1(x):
    return np.power(x, 2)

def f2(x):
    return np.abs(x)

def f3(x):
    
    if type(x) == np.ndarray:        
        result = []
        for num in x:
            result.append(f3(num))
        return result
    else:
        if (x < 1.0/3.0):
            return 1/2 * np.exp(x)
        else:
            return np.exp(x)


h_list = []
a = -1
b = 1

real_value_i1 = integrate.quad(f1, a, b)
real_value_i2 = integrate.quad(f2, a, b)
real_value_i3 = integrate.quad(f3, a, b)

error_1 = []
error_2 = []
error_3 = []

for i in range(1,21):
    h_list.append(np.exp2(-i))

for h in h_list:

    knotes = []
    h_ = a + h
    while( h_ < b):
        knotes.append(h_)
        h_ = h_+h
    
    values_1 = f1(np.asarray(knotes))
    values_2 = f2(np.asarray(knotes))
    values_3 = f3(np.asarray(knotes))
    #print(sum(values_2))
    trapez_quad_1 = h * (1.0/2.0*f1(a) + sum(values_1) + 1.0/2.0*f1(b))
    trapez_quad_2 = h * (1.0/2.0*f2(a) + sum(values_2) + 1.0/2.0*f2(b))
    trapez_quad_3 = h * (1.0/2.0*f3(a) + sum(values_3) + 1.0/2.0*f3(b))
    #print(trapez_quad_2)
    error_1.append(abs(trapez_quad_1 - real_value_i1[0]))
    error_2.append(abs(trapez_quad_2 - real_value_i2[0]))
    error_3.append(abs(trapez_quad_3 - real_value_i3[0]))

plt.loglog(h_list, error_1)
plt.loglog(h_list, error_2)
plt.loglog(h_list, error_3)
plt.show()
