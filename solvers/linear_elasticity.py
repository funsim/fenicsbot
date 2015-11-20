"""
linear elasticity
"""

from dolfin import *
from base_solver import BaseSolver

class LinearElasticitySolver(BaseSolver):

    @staticmethod
    def default_parameters():
        return {"domain": "UnitSquare",  # dimension
                "f": None,   # forcing term
                "E":10,
                "nu":0.3
                }

    # @staticmethod
    # def parameter_parsers():
    #     return {"domain": [], #no conversion - will default to lambda s: s
    #             "f": [lambda f: Constant(f.split(",")),
    #                   lambda f: Expression(f.split(","))],
    #             "E": [lambda E: Constant(E),
    #                   lambda E: Expression(E)], 
    #             "nu": [lambda nu: Constant(nu),
    #                   lambda nu: Expression(nu)]                                              
    #     }

    def solve(self):        
        f = self.params["f"]
        mesh = self.get_mesh()

        # Now that we know the dimension we can apply a default forcing if not
        # specified
        if f is None:
            n = mesh.geometry().dim()
            f = ("1,"*n)[:-1]

        # f = self.s2d(f)
        # E = self.s2d(self.params["E"])
        # nu = self.s2d(self.params["nu"])
        f, E, nu = map(self.s2d, 
                       [f, self.params["E"], self.params["nu"]])


        mu = E / (2*(1 + nu))
        lmbda = E*nu / ((1 + nu)*(1 - 2*nu))


        V = VectorFunctionSpace(mesh, 'CG', 1)

        # Define boundary conditions
        def u0_boundary(x, on_boundary):
            return on_boundary

        zero = Constant([0]*mesh.geometry().dim())
        bc = DirichletBC(V, zero, u0_boundary)

        # Define variational problem
        u = TrialFunction(V)
        v = TestFunction(V)

        def epsilon(v):
            return 0.5*(grad(v) + (grad(v)).T)
        def sigma(v):
            return 2*mu*epsilon(v) + lmbda*tr(epsilon(v))*Identity(len(v))

        a = inner(grad(v), sigma(u))*dx
        L = dot(v, f)*dx

        # Compute solution
        u = Function(V)
        solve(a == L, u, bc)
        
        self.solution = u

if __name__ == "__main__":

    params = LinearElasticitySolver.default_parameters()
    params["domain"] = "UnitSquare"
    # params["f"] = "1,sin(x[1])*x[0]"
    params["f"] = None
    params["E"] = 10
    params["nu"] = 0.3
    solver = LinearElasticitySolver(params)
    solver.solve()
    solver.plot()
