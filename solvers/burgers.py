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
                "nu": "0.00001",     # viscosity
                "domain": "UnitInterval",
        }

    def solve(self):

        mesh = self.get_mesh()
        self.update_parameters(self.params)

        T = float(self.params["dt"])
        timestep = self.s2d(self.params["dt"])
        f = self.s2d(self.params["f"])
        ic = self.s2d(self.params["ic"])
        nu = self.s2d(self.params["nu"])

        V = FunctionSpace(mesh, 'Lagrange', 1)

        def Dt(u, u_, timestep):
            return (u - u_)/timestep

        u_ = Function(ic)
        u = Function(V)
        v = TestFunction(V)

        timestep = Constant(dt)

        F = (Dt(u, u_, timestep)*v
             + u*u.dx(0)*v + nu*u.dx(0)*v.dx(0))*dx
        bc = DirichletBC(V, 0.0, "on_boundary")

        t = 0.0
        T = 0.2
        while (t <= T):
            solve(F == 0, u, bc)
            u_.assign(u)

            t += float(timestep)

        self.solution = u_

if __name__ == "__main__":

    params = BurgersSolver.default_parameters()

    solver = BurgersSolver(params)
    solver.solve()
    print solver.plot()
