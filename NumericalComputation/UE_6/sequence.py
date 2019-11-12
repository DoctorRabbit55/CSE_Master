import numpy as np
import math
import matplotlib.pyplot as plt


lim = 30
u = np.zeros(lim+1)
u[1] = 2.0  

for k in range(2,lim+1):    
    u[k] = 2**(k-1) * math.sqrt(2*(1-math.sqrt(1-(1/np.power(2,(k-1))*u[k-1])**2)))

        
error = abs(math.pi-u[1:])
plt.semilogy(np.arange(1,lim+1),error)
plt.show()
