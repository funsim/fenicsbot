"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Simplest example of computation and visualization with FEniCS.

-Laplace(u) = f on the unit square.
u = u0 on the boundary.
u0 = u = 1 + x^2 + 2y^2, f = -6.
"""

from base_solver import BaseSolver
from dolfin import *


class PoissonSolver(BaseSolver):

    @staticmethod
    def default_parameters():
        return {"f": "0",   # forcing term
                "domain": "UnitSquare", 

                ##BCs
                "bdy00": "0",
                "bdy01": "0",
                "bdy10": "0",
                "bdy11": "0",
                "bdy20": "0",
                "bdy21": "0"
        }

    def solve(self):
        mesh, bdys, facet_func = self.get_mesh(return_bdys=True)
        self.update_parameters(self.params)

        f = self.s2d(self.params["f"])

        V = FunctionSpace(mesh, 'Lagrange', 1)
        
        # Define boundary conditions
        #u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')
        # u_bc = Constant(0.0)
        # def u0_boundary(x, on_boundary):
        #         return on_boundary

        # bc = DirichletBC(V, u_bc, u0_boundary)
        bdy_names = self.subdomain_ordering()
        bcs = []
        for k in range(len(bdys)):
            bc_expr = self.s2d(self.params[bdy_names[k]])
            bcs.append(DirichletBC(V, bc_expr, facet_func, k))
            
        # Define variational problem
        u = TrialFunction(V)
        v = TestFunction(V)
        a = inner(nabla_grad(u), nabla_grad(v))*dx
        L = f*v*dx

        # Compute solution
        u = Function(V)
        solve(a == L, u, bcs)
        self.solution = u
        # print u.vector().array()

if __name__ == "__main__":

    params = PoissonSolver.default_parameters()
    params["domain"] = "UnitSquare"
    params["f"] = "0"
    params["bdy00"]="0"
    params["bdy01"]="1"

    solver = PoissonSolver(params)
    solver.solve()
    print solver.plot()
