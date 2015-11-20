from dolfin import *
from base_solver import BaseSolver

class StokesSolver(BaseSolver):

    @staticmethod
    def default_parameters():
        return {"domain": 2,
                "f": 0}

    def __init__(self, params):
        self.params = params
        self.update_parameters(params)

    def update_parameters(self, new_params):
        if "domain" in new_params:
            self.params["domain"] = str(new_params["domain"])
        if "f" in new_params:
            val = new_params["f"].split(",")
            self.params["f"] = Constant(val)

    def solve(self):
        f = self.params["f"]

        mesh = self.get_mesh()

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
    params["domain"] = "UnitCube"
    params["f"] = "2,3,4" #Expression(("x[0]*x[0]", "x[1]*x[1]"))

    solver = StokesSolver(params)
    solver.solve()
    solver.plot()
