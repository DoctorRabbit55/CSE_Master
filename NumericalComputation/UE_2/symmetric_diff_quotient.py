import numpy as np
from matplotlib import pyplot as plt

def neville_scheme(knotes, f, x):

    n = len(knotes)

    p = np.zeros((n,n))

    for k in range(n):
        for i in range(n-k):
            if k == 0:
                p[k, i] = f[i]
            else:
                p[k, i] = ((x - knotes[i+k])*p[k-1, i] - (x - knotes[i])*p[k-1, i+1])/(knotes[i]-knotes[i+k])

    return np.transpose(p)

def f(x):
    return np.tan(x)

def D_sym(h, x_0):
    h = np.asarray(h)
    return (f(x_0 + h) - f(x_0 - h)) / (h*2)

def plotConvergenceRates(knotes, p, n, true_value):
    plt.loglog(knotes[0:n], abs(p[0:n,0]-true_value), label="column 1")
    plt.loglog(knotes[0:n-1], abs(p[0:n-1,1]-true_value), label="column 2")
    plt.loglog(knotes[0:n-2], abs(p[0:n-2,2]-true_value), label="column 3")

    plt.loglog(knotes[0:n], np.power(knotes, 2), linestyle=":", color='cornflowerblue')
    plt.loglog(knotes[0:n], np.power(knotes, 3), linestyle=":", color='cornflowerblue')
    plt.loglog(knotes[0:n], np.power(knotes, 4), linestyle=":", color='cornflowerblue')

    plt.xlabel('min(knotes)')
    plt.ylabel('absolute error')
    plt.legend(loc="upper left")

    plt.show()

n = 8
x_0 = 0
knotes = []

for i in range(n):
    knotes.append(np.exp2(-i))

values = D_sym(knotes, x_0*np.ones(n))

p = neville_scheme(knotes, values, x_0)
print(p)
input('Press enter to continue..')

plotConvergenceRates(knotes, p, n, 1)


knotes_2 = np.power(knotes, 2)
p = neville_scheme(knotes_2, values, x_0)

print(p)
input('Press enter to continue..')

plotConvergenceRates(knotes, p, n, 1)

#-------------------------------------------------------------------

def f(x):
    x = np.asarray(x)
    print(x)
    print(np.power(np.amax([x, 0*x], axis=0), 3.0/2.0))
    return np.power(np.amax([x, 0*x], axis=0), 3.0/2.0)

values = D_sym(knotes, x_0*np.ones(n))

p = neville_scheme(knotes, values, x_0)
print(p)
input('Press enter to continue..')

plotConvergenceRates(knotes, p, n, 0)


knotes_2 = np.power(knotes, 2)
p = neville_scheme(knotes_2, values, x_0)

print(p)
input('Press enter to continue..')

plotConvergenceRates(knotes, p, n, 0)
