from base_solver import BaseSolver
from dolfin import *

class BurgersSolver(BaseSolver):

    @staticmethod
    def default_parameters():
        return {
                "f": "0",            # forcing term
                "ic": "sin(x[0])",   # initial condition
                "dt": "0.1",         # timestep
                "T": "1.0",          # final time
                "nu": "0.01",     # viscosity
                "domain": "UnitInterval",
        }

    def solve(self):

        mesh = self.get_mesh()
        self.update_parameters(self.params)

        T = float(self.params["T"])
        timestep = self.s2d(self.params["dt"])
        f = self.s2d(self.params["f"])
        ic = self.s2d(self.params["ic"])
        nu = self.s2d(self.params["nu"])

        if float(T)/float(timestep) > 50:
            raise ValueError, "This problem is a little too large for me."

        V = FunctionSpace(mesh, 'Lagrange', 1)

        def Dt(u, u_, timestep):
            return (u - u_)/timestep

        u_ = interpolate(ic, V)
        u = Function(V)
        v = TestFunction(V)

        F = (Dt(u, u_, timestep)*v
             + u*u.dx(0)*v + nu*u.dx(0)*v.dx(0))*dx
        # bc = DirichletBC(V, 0.0, "on_boundary")
        bcs = self.get_bcs(V)

        t = 0.0
        while (t <= T):
            solve(F == 0, u, bcs)
            u_.assign(u)

            t += float(timestep)

        self.solution = u_

if __name__ == "__main__":

    params = BurgersSolver.default_parameters()

    solver = BurgersSolver(params)
    solver.solve()
    print solver.plot()
