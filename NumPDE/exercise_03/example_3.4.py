from ngsolve import *
import netgen.geom2d as geom2d
from netgen.geom2d import unit_square

import numpy as np
from matplotlib import pyplot as plt

ngsglobals.msg_level = 1

geo = geom2d.SplineGeometry()

p1, p2, p3, p4 = [ geo.AppendPoint(x,y) for x,y in [(0.0,0.0), (1.0,0.0), (1.0, 1.0), (0.0,1.0)]]
p5, p6, p7, p8 = [ geo.AppendPoint(x,y) for x,y in [(0.3,0.5), (0.5,0.5), (0.5, 0.7), (0.3,0.7)]]

geo.Append(["line", p1, p2], leftdomain=2, rightdomain=0, bc="1")
geo.Append(["line", p2, p3], leftdomain=2, rightdomain=0, bc="2")
geo.Append(["line", p3, p4], leftdomain=2, rightdomain=0, bc="3")
geo.Append(["line", p4, p1], leftdomain=2, rightdomain=0, bc="4")

geo.Append(["line", p5, p6], leftdomain=1, rightdomain=2, bc="5")
geo.Append(["line", p6, p7], leftdomain=1, rightdomain=2, bc="6")
geo.Append(["line", p7, p8], leftdomain=1, rightdomain=2, bc="7")
geo.Append(["line", p8, p5], leftdomain=1, rightdomain=2, bc="8")

geo.SetMaterial (2, "outer")
geo.SetMaterial (1, "inner")

mesh = Mesh(geo.GenerateMesh(maxh=0.05))

fes = H1(mesh, order=2, dirichlet=[1,2,3,4])

u = fes.TrialFunction()
v = fes.TestFunction()

f_co = CoefficientFunction([1,0])
f = LinearForm(fes)
f += f_co*v*dx

lam = CoefficientFunction([10,1])
a = BilinearForm(fes, symmetric=False)
a += lam*grad(u)*grad(v)*dx

a.Assemble()
f.Assemble()

gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

Draw(gfu)
Draw(lam*grad(gfu), mesh, "Flux")

space_flux = HDiv(mesh, order=2)
gf_flux = GridFunction(space_flux)

flux = lam*grad(gfu)
gf_flux.Set(flux)

# dx
#help(Integrate)
flux_y = Integrate(gf_flux[1], mesh, definedon=mesh.Boundaries("5|7"), region_wise=True)
flux_x = Integrate(gf_flux[0], mesh, definedon=mesh.Boundaries("6|8"), region_wise=True)

print("x: ",flux_x)
print("y: ",flux_y)

print("flux: ", flux_x[5]-flux_x[7]+flux_y[6]-flux_y[4])


