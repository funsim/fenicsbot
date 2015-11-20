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


class PoissonSolver(BaseSolver):

    @staticmethod
    def default_parameters():
        return {"f": Constant(0),   # forcing term
                "domain": "UnitSquare" #
                }

    def __init__(self, params):
        self.params = PoissonSolver.default_parameters()
        self.update_parameters(params)

    def update_parameters(self, new_params):
        if "domain" in new_params:
            self.params["domain"] = str(new_params["domain"])
        if "f" in new_params:
            self.params["f"] = Constant(float(new_params["f"]))

    def solve(self):
        f = self.params["f"]

        mesh = self.get_mesh()

        V = FunctionSpace(mesh, 'Lagrange', 1)

        # Define boundary conditions
        #u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')
        u0 = Constant(0.0)
        def u0_boundary(x, on_boundary):
                return on_boundary

        bc = DirichletBC(V, u0, u0_boundary)

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
    params["f"] = 2 # Expression("x[0]*x[0]")

    solver = PoissonSolver(params)
    solver.solve()
    print solver.plot()
