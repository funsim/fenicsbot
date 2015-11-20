from dolfin import *
from base_solver import BaseSolver

class StokesSolver(BaseSolver):

    @staticmethod
    def default_parameters():
        return {"domain": "UnitSquare",
                "f": None
               }

    def __init__(self, params):
        self.params = StokesSolver.default_parameters()
        self.update_parameters(params)

    def update_parameters(self, new_params):
        if "domain" in new_params:
            self.params["domain"] = str(new_params["domain"])
        if "f" in new_params:
            self.params["f"] = new_params["f"]

    def solve(self):
        f = self.params["f"]
        mesh = self.get_mesh()

        # Now that we know the dimension we can apply a default forcing if not
        # specified
        if f is None:
            f = ",".join(["1"]*mesh.geometry().dim())

        fvals = f.split(",")
        try:
            f = Constant(fvals)
        except:
            f = Expression(fvals)

        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)
        W = V * Q

        # Define boundary conditions
        #u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')
        def dirichlet_boundary(x, on_boundary):
            return x[0] > 1.0 - DOLFIN_EPS or x[0] < DOLFIN_EPS

        zero = Constant([0]*mesh.geometry().dim())
        bcs = DirichletBC(W.sub(0), zero, dirichlet_boundary)

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

        self.solution = u

if __name__ == "__main__":

    params = StokesSolver.default_parameters()
    params["domain"] = "UnitSquare"
    params["f"] = "1,sin(x[1])*x[0]"

    solver = StokesSolver(params)
    solver.solve()
    solver.plot()
