import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def f1(x):
    return np.power(x, 0.2)

def f2(x):
    return np.power(x, 10)

def f3(x):  
    return np.power(x, 2)

def neville_scheme(knotes, f, x):

    n = len(knotes)

    p = np.zeros((n,n))

    for k in range(n):
        for i in range(n-k):
            if k == 0:
                p[k, i] = f[i]
            else:
                p[k, i] = ((x - knotes[i+k])*p[k-1, i] - (x - knotes[i])*p[k-1, i+1])/(knotes[i]-knotes[i+k])

    return p


h_list = []
a = 0
b = 1

real_value_1 = integrate.quad(f1, a, b)[0]
real_value_2 = integrate.quad(f2, a, b)[0]
real_value_3 = integrate.quad(f3, a, b)[0]

trapez_quad_1 = []
trapez_quad_2 = []
trapez_quad_3 = []

I = 5

for i in range(1,I):
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
    
    trapez_quad_1.append(h * (1.0/2.0*f1(a) + sum(values_1) + 1.0/2.0*f1(b)))
    trapez_quad_2.append(h * (1.0/2.0*f2(a) + sum(values_2) + 1.0/2.0*f2(b)))
    trapez_quad_3.append(h * (1.0/2.0*f3(a) + sum(values_3) + 1.0/2.0*f3(b)))


table_1 = neville_scheme(np.power(h_list, 2), trapez_quad_1, 0)
table_2 = neville_scheme(np.power(h_list, 2), trapez_quad_2, 0)
table_3 = neville_scheme(np.power(h_list, 2), trapez_quad_3, 0)

plt.loglog(h_list[0:I-1], abs(table_1[0:I-1,0]-np.ones(I-1)*real_value_1), label="column 1")
plt.loglog(h_list[0:I-2], abs(table_1[0:I-2,1]-np.ones(I-2)*real_value_1), label="column 2")
plt.loglog(h_list[0:I-3], abs(table_1[0:I-3,2]-np.ones(I-3)*real_value_1), label="column 3")
plt.loglog(h_list, np.power(h_list, 1), linestyle=":", color='cornflowerblue', label='h^1')
plt.loglog(h_list, np.power(h_list, 4), linestyle=":", color='cornflowerblue', label='h^4')

plt.legend()
plt.show()

plt.loglog(h_list[0:I-1], abs(table_2[0:I-1,0]-np.ones(I-1)*real_value_2), label="column 1")
plt.loglog(h_list[0:I-2], abs(table_2[0:I-2,1]-np.ones(I-2)*real_value_2), label="column 2")
plt.loglog(h_list[0:I-3], abs(table_2[0:I-3,2]-np.ones(I-3)*real_value_2), label="column 3")
plt.loglog(h_list, np.power(h_list, 1), linestyle=":", color='cornflowerblue', label='h^1')
plt.loglog(h_list, np.power(h_list, 4), linestyle=":", color='cornflowerblue', label='h^4')

plt.legend()
plt.show()

plt.loglog(h_list[0:I-1], abs(table_3[0:I-1,0]-np.ones(I-1)*real_value_3), label="column 1")
plt.loglog(h_list[0:I-2], abs(table_3[0:I-2,1]-np.ones(I-2)*real_value_3), label="column 2")
plt.loglog(h_list[0:I-3], abs(table_3[0:I-3,2]-np.ones(I-3)*real_value_3), label="column 3")
plt.loglog(h_list, np.power(h_list, 8), linestyle=":", color='cornflowerblue', label='h^8')
plt.loglog(h_list, np.power(h_list, 12), linestyle=":", color='cornflowerblue', label='h^12')

plt.legend()
plt.show()
