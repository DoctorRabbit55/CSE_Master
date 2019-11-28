import numpy as np
import numpy.linalg as linalg
import random

N = np.arange(2,11)

for i in range(len(N)):

    A = np.zeros((N[i],2))
    y = np.zeros((N[i],1))

    for idx in range(N[i]):
        x = random.random()*1/(float)(len(N)) - 1/(float)(len(N))
        
        A[idx, 0] = x
        A[idx, 1] = x**3

        y[idx] = np.sin(x)

    X = np.dot(A.T,A)
    b = np.dot(A.T,y)

    sol = linalg.solve(X, b)
       
    print(sol)


