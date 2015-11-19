"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Simplest example of computation and visualization with FEniCS.

-Laplace(u) = f on the unit square.
u = u0 on the boundary.
u0 = u = 1 + x^2 + 2y^2, f = -6.
"""

from dolfin import *

#take as input the type of domain


class StokesSolver(object):

    @staticmethod
    def default_parameters():
        return {"D": 2,
                "f": 0}

    def __init__(self, params):
        self.params = params

    def __call__(self):
        D = params["D"]
        f = params["f"]

        if D == 1:
            mesh = UnitIntervalMesh(20)
            zero = Constant(0.0)
        elif D==2:
            mesh = UnitSquareMesh(20, 20)
            zero = Constant((0.0,0.0))
        elif D==3:
            mesh = UnitCubeMesh(20, 20, 20)
            zero = Constant((0.0,0.0,0.0))

        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)
        W = V * Q

        # Define boundary conditions
        #u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')
        def Dirichlet_boundary(x, on_boundary):
            return x[1] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS

        bcs = DirichletBC(W.sub(0), zero, Dirichlet_boundary)

        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        a = inner(grad(u), grad(v))*dx + div(v)*p*dx + q*div(u)*dx
        L = inner(f, v)*dx

        # Form for use in constructing preconditioner matrix
        b = inner(grad(u), grad(v))*dx + p*q*dx

        # Assemble system
        A, bb = assemble_system(a, L, bcs)

        # Assemble preconditioner system
        P, btmp = assemble_system(b, L, bcs)

        # Create Krylov solver and AMG preconditioner
        solver = KrylovSolver("tfqmr", "amg")

        # Associate operator (A) and preconditioner matrix (P)
        solver.set_operators(A, P)

        # Solve
        U = Function(W)
        solver.solve(U.vector(), bb)

        # Get sub-functions
        u, p = U.split()

        # Save solution in VTK format
        solution_u = plot(u)
        solution_u.write_png("solution_u")

        solution_p = plot(p)
        solution_p.write_png("solution_p")

        plot_mesh = plot(mesh)
        plot_mesh.write_png("mesh")

if __name__ == "__main__":

    params = StokesSolver.default_parameters()
    params["D"] = 2
    params["f"] = Expression(("x[0]*x[0]","x[1]*x[1]"))

    solver = StokesSolver(params)
    solver()
