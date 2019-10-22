import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def f(x):
    return 1/(4-np.power(x,2))

def neville_scheme(knotes, values):

    n = len(knotes)
    p = sp.zeros(n)

    for k in range(n):
        for i in range(n-k):
            if k == 0:
                p[k, i] = values[i]
            else:
                p[k, i] = sp.simplify(((x - knotes[i+k])*p[k-1, i] - (x - knotes[i])*p[k-1, i+1])/(knotes[i]-knotes[i+k]))

    return sp.transpose(p)

max_error_list = [];

for n in range(1,10):

    a = -1
    b = 1
    x = sp.symbols('x')     
    knotes = []
    values = []

    for i in range(n):
        knotes.append((a+b)/2 + (b-a)/2 * np.cos(np.pi*(2*i+1)/(2*(n-1)+2)))

    values = f(knotes)
    
    p = neville_scheme(knotes, values)
    p = sp.simplify(p[0,n-1])

    uni_points = np.arange(-1, 1, 2/100)
    app_values = []
    real_values = []
    for i in range(100):
        app_values.append(p.subs(x, uni_points[i]))
        real_values.append(f(uni_points[i]))

    max_error_list.append(np.amax(abs(np.asarray(real_values) - np.array(app_values))))

    #plt.plot(uni_points, app_values)
    #plt.plot(uni_points, real_values)
    #plt.show()

plt.semilogy(np.arange(1,10), max_error_list)
plt.show()
