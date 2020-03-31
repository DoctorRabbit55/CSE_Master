from ngsolve import *
import netgen.geom2d as geom2d

#.. define mesh...
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

geo.SetMaterial (1, "inner")
geo.SetMaterial (2, "outer")

mesh = Mesh(geo.GenerateMesh(maxh=0.1))

order_u = 2

# define discrete conforming spaces of Sigma and V
# any ideal why we need a different order here?
# try to use different combinations later as well
Sigma = HDiv(mesh, order = order_u + 1)
V = L2(mesh, order = order_u, dirichlet=[1,2,3,4])

#This defines a "compound FESpace"
X = FESpace([Sigma, V])

#This gives you a trial function sigma in Sigma and u in V
(sigma, u) = X.TrialFunction()
(tau, v) = X.TestFunction()
#...

# define a "big" bilinearform on X which includes all the integrals of (4)
# you can use div(sigma) to get the divergence
lam = CoefficientFunction([1,10])
B = BilinearForm(X)
B += div(sigma)*v*dx
B += -1/lam*sigma*tau*dx + div(tau)*u*dx
#...

# and a linearform on X
f_co = CoefficientFunction([1,0])
F = LinearForm(X)
F += f_co*v*dx
#...

B.Assemble()
F.Assemble()

# solution on X
sol = GridFunction(X)
sigma_sol = sol.components[0] #this gives you the sigma solution
u_sol = sol.components[1] #this gives you the u solution

# solve the problem
# NOTE: the Bilinearform B is NOT coercive -> not SPD, hence we can not use a sparsecholesky solver
# use: ... B.mat.Inverse(X.FreeDofs(), inverse = "umfpack")
sol.vec.data = B.mat.Inverse(X.FreeDofs(), inverse = "umfpack") * F.vec

Draw(sol.components[0], mesh, "sol")

#normal vector needed for the fluxes
n = specialcf.normal(2)

# evaluate the flux using the sigma solution!
# Similar as for the H1, the var framework from the cont setting transfers to the discrete setting
# For sigma you can now evaluate sigma * n (since the normal trace is a cont operator for the H(div)

flux = -sol.components[0] * n
print("flux = ", Integrate(flux, mesh, definedon=mesh.Boundaries("5|6|7|8")))
