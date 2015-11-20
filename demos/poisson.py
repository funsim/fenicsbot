"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Simplest example of computation and visualization with FEniCS.

-Laplace(u) = f on the unit square.
u = u0 on the boundary.
u0 = u = 1 + x^2 + 2y^2, f = -6.
"""

import tempfile
from base_solver import BaseSolver
from dolfin import *
# from base_solver.BaseSolver import string_to_dolfin as s2d


class PoissonSolver(BaseSolver):

    @staticmethod
    def default_parameters():
        return {"f": "0",   # forcing term
                "domain": "UnitSquare" #
        }

    def solve(self):
        mesh = self.get_mesh()
        self.update_parameters(self.params)

        f = self.s2d(self.params["f"])

        V = FunctionSpace(mesh, 'Lagrange', 1)

        # Define boundary conditions
        #u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')
        u_bc = Constant(0.0)
        def u0_boundary(x, on_boundary):
                return on_boundary

        bc = DirichletBC(V, u_bc, u0_boundary)

        # Define variational problem
        u = TrialFunction(V)
        v = TestFunction(V)
        a = inner(nabla_grad(u), nabla_grad(v))*dx
        L = f*v*dx

        # Compute solution
        u = Function(V)
        solve(a == L, u, bc)
        self.solution = u

if __name__ == "__main__":

    params = PoissonSolver.default_parameters()
    params["domain"] = "Dolfin"
    params["f"] = "x[0]*x[1]"

    solver = PoissonSolver(params)
    solver.solve()
    print solver.plot()
