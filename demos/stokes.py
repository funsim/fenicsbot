"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Simplest example of computation and visualization with FEniCS.

-Laplace(u) = f on the unit square.
u = u0 on the boundary.
u0 = u = 1 + x^2 + 2y^2, f = -6.
"""

from dolfin import *

#take as input the type of domain


def stokes_solver(D, f):

	if D == 1:
		mesh = UnitIntervalMesh(20)
	elif D==2:
	   	mesh = UnitSquareMesh(20, 20)
	elif D==3:
		mesh = UnitCubeMesh(20, 20, 20)

	V = VectorFunctionSpace(mesh, "CG", 2)
	Q = FunctionSpace(mesh, "CG", 1)
	W = V * Q

	# Define boundary conditions
	#u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')

	u0 = Constant(0.0)
	def u0_boundary(x, on_boundary):
		return on_boundary &&  x[1] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS

	bc = DirichletBC(V, u0, u0_boundary)

	# Define variational problem
	u = TrialFunction(V)
	v = TestFunction(V)
#	f = Constant(-6.0)
	a = inner(nabla_grad(u), nabla_grad(v))*dx
	L = f*v*dx

	# Compute solution
	u = Function(V)
	solve(a == L, u, bc)

	# Plot solution and mesh
	solution = plot(u)
	solution.write_png("solution")
	plot(mesh)

	plot_mesh = plot(mesh)
	plot_mesh.write_png("mesh")

if __name__ == "__main__":
    D = 2
    f = Expression("x[0]*x[0]")
    poisson_solver(D, f)
