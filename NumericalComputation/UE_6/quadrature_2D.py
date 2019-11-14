import numpy as np
from scipy import integrate

import matplotlib.pyplot as plt



def Q(a,b,c,d,f):  
    return (b-a)*(d-c)*f((a+b)/2,(c+d)/2)


def midpoint(a,b,c,d,f,tau,hmin):
    
    if((b-a) <= hmin or (d-c) <= hmin):
        return Q(a,b,c,d,f)
    
    finerScale = Q(a,(a+b)/2,c,(c+d)/2,f) + Q((a+b)/2,b,c,(c+d)/2,f) + Q((a+b)/2,b,(c+d)/2,d,f) + Q(a,(a+b)/2,(c+d)/2,d,f)
    
    if abs(Q(a,b,c,d,f) - finerScale) <= tau:
        return Q(a,b,c,d,f)
    else:
        return midpoint(a,(a+b)/2,c,(c+d)/2,f,tau/2,hmin) + midpoint((a+b)/2,b,c,(c+d)/2,f,tau/2,hmin) + midpoint((a+b)/2,b,(c+d)/2,d,f,tau/2,hmin) + midpoint(a,(a+b)/2,(c+d)/2,d,f,tau/2,hmin)


def f1(x,y):  
    return x**2


def f2(x,y):
    
    if(x<y):
        return 0
    else:
        return 1
    
    
        
i = np.arange(0,16)
a = 0
b = 1

real_value, error = integrate.dblquad(f1, a, b, a, b)
error_list = np.zeros(16)
print(real_value)
tau = 1/np.power(2,i)

for j in i:
    error_list[j] = abs(midpoint(a, b, a, b, f1, tau[j], 0.1) - real_value)

plt.loglog(tau, error_list)
plt.show()

#real_value, error = integrate.dblquad(f2, a, b, a, b)
real_value = 1/2
print(real_value)

for j in i:
    error_list[j] = abs(midpoint(a, b, a, b, f2, tau[j], 0.001) - real_value)
    
plt.loglog(tau, error_list)
plt.show()
    
