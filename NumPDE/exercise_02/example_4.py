from ngsolve import *
import netgen.geom2d as geom2d
from netgen.geom2d import unit_square

import numpy as np
from matplotlib import pyplot as plt

ngsglobals.msg_level = 1

k = 3 # change this parameter for solutions of different order

# CONCLUSION
# smaller mesh will not solve problem accuratly -> higher order solution is required

exact_value_left = 1/6
exact_value_right = -1/3

error_list = []

h_list = []
for i in range(1,8):
  h_list.append(2**(-i))

for h in h_list:

  mesh = Mesh(unit_square.GenerateMesh(maxh=h))

  fes = H1(mesh, order=k, dirichlet=[2,4])

  u = fes.TrialFunction()
  v = fes.TestFunction()

  f = LinearForm(fes)
  f += x * v * dx

  a = BilinearForm(fes, symmetric=True)
  a += grad(u)*grad(v)*dx

  a.Assemble()
  f.Assemble()

  gfu = GridFunction(fes)
  gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

  dudx = GridFunction(fes)
  dudx.Set(grad(gfu)[0])

  points = np.arange(0, 1+h, h)

  integral_left = 0
  integral_right = 0

  for i in range(1,len(points)):
    point_left_1 = mesh(0, points[i-1])
    point_left_2 = mesh(0, points[i])

    point_right_1 = mesh(1, points[i-1])
    point_right_2 = mesh(1, points[i])

    integral_left += h/2*(dudx(point_left_1) + dudx(point_left_2))
    integral_right += h/2*(dudx(point_right_1) + dudx(point_right_2))

  print("numeric results for h=", h)
  print("W_l:  ", integral_left)
  print("W_r: ", integral_right)
  print("W_l + W_r: ", integral_left + integral_right)
  print("---------------")

  integrals = Integrate(dudx, mesh, BND, region_wise=True)
  print("intern function for h=", h)
  print("W_l: ", integrals[3])
  print("W_r: ", integrals[1])
  print("W_l + W_r: ", integrals[3] + integrals[1])

  error_list.append((abs(integral_left - exact_value_left), abs(integral_right - exact_value_right) \
  , abs(integral_right - exact_value_right) + abs(integral_right - exact_value_right)))


plt.loglog(h_list, [x[0] for x in error_list], label='left')
plt.loglog(h_list, [x[1] for x in error_list], label='right')
plt.loglog(h_list, [x[2] for x in error_list], label='added')
plt.legend()
plt.show()
