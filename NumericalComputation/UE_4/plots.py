import numpy as np
from matplotlib import pyplot as plt

C = 2
alpha = 1/2

h = np.arange(1,100)

plt.loglog(h, C*np.power(h,alpha))
plt.show()

plt.semilogy(h, C*np.exp(-alpha*h))
plt.show()

plt.semilogx(h, C*np.exp(-alpha/h))
plt.show()
