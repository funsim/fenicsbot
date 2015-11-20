from dolfin import *
from base_solver import BaseSolver

class StokesSolver(BaseSolver):

    @staticmethod
    def default_parameters():
        return {"domain": "UnitSquare",
                "f": None
        }


    def solve(self):
        f = self.params["f"]
        mesh = self.get_mesh()

        # Now that we know the dimension we can apply a default forcing if not
        # specified
        if f is None:
            n = mesh.geometry().dim()
            f = ("1,"*n)[:-1]

        f = self.s2d(f)

        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)
        W = V * Q

        # Define boundary conditions
        #u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')
        # def dirichlet_boundary(x, on_boundary):
        #     return x[0] > 1.0 - DOLFIN_EPS or x[0] < DOLFIN_EPS

        # zero = Constant([0]*mesh.geometry().dim())
        zero = ("0,"*mesh.geometry().dim())[:-1]
        
        # bcs = DirichletBC(W.sub(0), zero, dirichlet_boundary)
        bcs = self.get_bcs(W.sub(0), default=zero)

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

        self.solution = interpolate(u, V)

if __name__ == "__main__":

    params = StokesSolver.default_parameters()
    #params["domain"] = "UnitSquare"
    # params["f"] = "1,sin(x[1])*x[0]"
    #params["f"] = None

    solver = StokesSolver(params)
    solver.solve()
    solver.plot()
