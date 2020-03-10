from ngsolve import *
import netgen.geom2d as geom2d
from netgen.geom2d import unit_square

from matplotlib import pyplot as plt

ngsglobals.msg_level = 0

geo = geom2d.SplineGeometry()
p1 = geo.AppendPoint(0,0)
p2 = geo.AppendPoint(2,0)
p3 = geo.AppendPoint(2,1)
p4 = geo.AppendPoint(1,1)
p5 = geo.AppendPoint(1,2)
p6 = geo.AppendPoint(0,2)

geo.Append(["line", p1, p2])
geo.Append(["line", p2, p3])
geo.Append(["line", p3, p4])
geo.Append(["line", p4, p5])
geo.Append(["line", p5, p6])
geo.Append(["line", p6, p1])

mesh = Mesh(geo.GenerateMesh(maxh=0.1))
#mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))

fes_exact = H1(mesh, order=10, dirichlet=[1,2,3,4,5,6])

u_exact = fes_exact.TrialFunction()
v_exact = fes_exact.TestFunction()

f_exact = LinearForm(fes_exact)
f_exact += v_exact * dx

a_exact = BilinearForm(fes_exact, symmetric=True)
a_exact += grad(u_exact)*grad(v_exact)*dx

a_exact.Assemble()
f_exact.Assemble()

gfu_exact = GridFunction(fes_exact)
gfu_exact.vec.data = a_exact.mat.Inverse(fes_exact.FreeDofs(), inverse="sparsecholesky") * f_exact.vec
grad_exact = grad(gfu_exact)

Draw(gfu_exact)
Draw(-grad(gfu_exact), mesh, "Flux")

data_list = []

for i in range(len(gfu_exact.vec.data)):
  data_list.append(gfu_exact.vec.data[i])

data_list.sort()
print("mean: ", data_list[int(len(data_list)/2)])


L2_error_list = []
H1_error_list = []

for k in range(1,9):

  fes = H1(mesh, order=k, dirichlet=[1,2,3,4,5,6])

  u = fes.TrialFunction()
  v = fes.TestFunction()

  f = LinearForm(fes)
  f += v * dx
  f.Assemble()

  a = BilinearForm(fes, symmetric=True)
  a += grad(u)*grad(v)*dx
  a.Assemble()

  gfu = GridFunction(fes)
  gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

  L2_error = sqrt(Integrate((gfu-gfu_exact)*(gfu-gfu_exact), mesh))
  L2_error_list.append(L2_error)
  print("L2_error: ", L2_error)

  grad_ = grad(gfu)
  H1_error = L2_error + sqrt(Integrate((grad_-grad_exact)*(grad_-grad_exact), mesh))
  H1_error_list.append(H1_error)
  print("H1_error: ", H1_error)

#Draw(gfu)
#raw(-grad(gfu), mesh, "Flux")

plt.plot(range(len(L2_error_list)), L2_error_list)
plt.show()


