import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt  

A0 = np.array([[3,2,1], [51,3,2], [3,0,0]])
lmax = 10
A_list = []
A_list.append(A0)

for i in range(lmax):
    Q,R = la.qr(A0)
    A0 = R*Q
    A_list.append(A0)

print(np.diagonal(A0))
