import numpy as np
import matplotlib.pyplot as plt

def f1(x):
    return np.power(x, 2)

def f2(x):
    return np.abs(x)

def f3(x):
    if (x < 1.0/3.0):
        return 1/2 * np.exp(x)
    else:
        return np.exp(x)


h_list = []
a = -1
b = 1
for i in range(1,21):
    h_list.append(np.exp2(-i))

for h in h_list:

    knotes = np.arange(a, b, h)
    
    values_1 = f1(knotes)
    values_2 = f2(knotes)
    values_3 = f3(knotes)

    trapez_quad_1 = h * (1.0/2.0*f1(a) + sum(values_1) + 1.0/2.0*f1(b))
    trapez_quad_2 = h * (1.0/2.0*f2(a) + sum(values_2) + 1.0/2.0*f2(b))
    trapez_quad_3 = h * (1.0/2.0*f3(a) + sum(values_3) + 1.0/2.0*f3(b))

    
