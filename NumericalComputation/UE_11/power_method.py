import numpy as np
import numpy.linalg as la

def power_method(A, x):
  x = A.dot(x)
  x = x / la.norm(x)
  ev = x.T.dot(A.dot(x))
  
  return x, ev
  
  
A = np.array([[2,0], [0,-2]])

x0 = np.array([1,0]).T
x0 = x0 / la.norm(x0)
ev = 0

for i in range(20):
  x0, ev = power_method(A, x0)
  
print("x: ", x0)
print("eigenvalue: ", ev)


