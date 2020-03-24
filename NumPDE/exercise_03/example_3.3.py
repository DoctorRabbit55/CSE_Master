from ngsolve import *
import netgen.geom2d as geom2d
from netgen.geom2d import unit_square

import numpy as np
from matplotlib import pyplot as plt

ngsglobals.msg_level = 1

mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))

fes = H1(mesh, order=3, dirichlet=[2,4])

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

# define test function
w_r = GridFunction(fes)
w_r.Set(x)

w_l = GridFunction(fes)
w_l.Set(1-x)

# Integrating
W_r = Integrate(grad(gfu)*grad(w_r)-x*w_r, mesh)
W_l = Integrate(grad(gfu)*grad(w_l)-x*w_l, mesh)

print("----------------------------------")
print("Integrating")
print("flux right: ", W_r)
print("flux left: ", W_l)

# Matrix multiplication
b = gfu.vec.CreateVector()
b.data = a.mat*gfu.vec - f.vec
W_r = InnerProduct(b, w_r.vec)
W_l = InnerProduct(b, w_l.vec)

print("----------------------------------")
print("Matrix multiplication")
print("flux right: ", W_r)
print("flux left: ", W_l)

