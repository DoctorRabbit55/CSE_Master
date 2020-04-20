from ngsolve import *
import netgen.geom2d as geom2d

#.. define mesh...
ngsglobals.msg_level = 1

geo = geom2d.unit_square

mesh = Mesh(geo.GenerateMesh(maxh=0.1))

order_u = 1
t = 0.00001

V1 = H1(mesh, order=order_u, dirichlet=[1,2,3,4])
V2 = H1(mesh, order=order_u, dirichlet=[1,2,3,4])
V3 = H1(mesh, order=order_u, dirichlet=[1,2,3,4])


X = FESpace([V1, V2, V3])


(betax, betay, w) = X.TrialFunction()
(deltax, deltay, v) = X.TestFunction()


B = BilinearForm(X)
B += (grad(betax)[0]*grad(deltax)[0]+grad(betay)[1]*grad(deltay)[1])*dx + \
     1/t**2*((grad(w)[0]-betax)*(grad(v)[0]-deltax)+(grad(w)[1]-betay)*(grad(v)[1]-deltay))*dx


F = LinearForm(X)
F += v*dx


B.Assemble()
F.Assemble()


sol = GridFunction(X)

sol.vec.data = B.mat.Inverse(X.FreeDofs(), inverse = "umfpack") * F.vec

Draw(sol.components[2], mesh, "w")
Draw(sol.components[0], mesh, "beta_x")
Draw(sol.components[1], mesh, "beta_y")

print("Integral[0] = ", Integrate(grad(sol.components[2])[0]-sol.components[0], mesh))
print("Integral[1] = ", Integrate(grad(sol.components[2])[1]-sol.components[1], mesh))
