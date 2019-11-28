import numpy as np
import numpy.linalg as linalg
import random

n = np.arange(2, 11)

for i in n:

    N = 2**i
    
    A = np.zeros((i,2))
    y = np.zeros((i,1))

    for idx in range(i):
        x = random.random() * 2/(float)(N) - 1/(float)(N)
        
        A[idx, 0] = x
        A[idx, 1] = x**3

        y[idx] = np.sin(x)

    X = np.dot(A.T,A)
    b = np.dot(A.T,y)

    sol = linalg.solve(X, b)
       
    print(sol)
