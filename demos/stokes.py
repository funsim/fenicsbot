from dolfin import *
from base_solver import BaseSolver

class StokesSolver(BaseSolver):

    @staticmethod
    def default_parameters():
        return {"D": 2,
                "f": 0}

    def __init__(self, params):
        self.params = params
        self.update_parameters(params)

    def update_parameters(self, new_params):
        if "D" in new_params:
            self.params["D"] = int(new_params["D"])
        if "f" in new_params:
            val = new_params["f"].split(",")
            self.params["f"] = Constant(val)

    def solve(self):
        D = self.params["D"]
        f = self.params["f"]

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

        self.solution = u

if __name__ == "__main__":

    params = StokesSolver.default_parameters()
    params["D"] = "2"
    params["f"] = "2, 4" #Expression(("x[0]*x[0]", "x[1]*x[1]"))

    solver = StokesSolver(params)
    solver.solve()
    solver.plot()
