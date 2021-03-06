import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import os
import tempfile
from dolfin import *
import matplotlib.pyplot as plt


class BaseSolver(object):

    @staticmethod
    def default_parameters():
        raise NotImplementedError

    def s2d(self, s):
        """ Attempt to convert a string to a dolfin Constant or Expression. """
        s = str(s)  # Need to convert to str because the Twitter module gives
                    # unicode strings which dolfin cannot handle properly
        if "," in s:
            s = s.split(",")

        try:
            return Constant(s)
        except:
            try:
                return Expression(s)
            except RuntimeError:

                # Support expressing coordaintes with x, y and z
                if isinstance(s, list):
                    s = [ss.replace("x", "x[0]").replace("y", "x[1]").replace("z",
                        "x[2]") for ss in s]
                else:
                    s = s.replace("x", "x[0]").replace("y", "x[1]").replace("z",
                        "x[2]")
                return Expression(s)


    def parameter_parsers(self):
        """
        Should be overridden by solver subclass, and return a dictionary of
        possible conversions in the format

        { "parname": [lambda s: conversion_1(s),
                      lambda s: conversion_2(s), ...]}

        a "default converter" lambda s: s is always used if no conversion works
        for example:

        {"f":    [lambda s: Constant(float(s)),
                  lambda s: Expression(f)],
         "mesh": []}

        would mean that parameters f are read as constants if possible,
        expressions if not (and as strings if that doesn't work either)
        while the mesh parameter is just stored as string without any conversion
        """

        ## specify arguments as a dict {argname: list of possible "conversions"}
        ## should be called by the solve() method, and not before

        raise NotImplementedError



    def solve(self):
        raise NotImplementedError

    def update_parameters(self, new_parameters):
        for parname in new_parameters:
            self.params[parname] = new_parameters[parname]

            # parsers = self.__class__.parameter_parsers()[parname]
            # parsers.append(lambda s: s) # default if nothing works
            # try:





    def __init__(self, params):
        self.params = self.__class__.default_parameters()
        self.update_parameters(params)
        self.canonicalize_params()

    def canonicalize_params(self):
        """ Apply some canonicalization the the parameters """
        domain = self.params["domain"].lower()

        if domain == "unitsquare" or domain == "square":
            self.params["domain"] = "UnitSquare"

        if domain == "unitinterval" or domain == "interval":
            self.params["domain"] = "UnitInterval"

        if domain == "unitcube" or domain == "cube":
            self.params["domain"] = "UnitCube"

        if domain == "unitcircle" or domain == "circle":
            self.params["domain"] = "UnitCircle"

        if domain == "l":
            self.params["domain"] = "L"

    def get_mesh(self):
        """
        Returns mesh given by self.params["domain"], and creates a
        division into subdomains as well as a corresponding
        FacetFunction. The latter two is stored as properties of self.
        """

        domain = self.params["domain"]
        here = os.path.dirname(__file__)

        if domain == "UnitInterval":
            mesh = UnitIntervalMesh(20)
        elif domain == "UnitSquare":
            mesh = UnitSquareMesh(20, 20)
        elif domain == "UnitCube":
            mesh = UnitCubeMesh(5, 5, 5)
        elif domain == "Dolfin":
            mesh = Mesh(os.path.join(here, "meshes", "dolfin.xml.gz"))
        elif domain == "Circle":
            mesh = Mesh(os.path.join(here, "meshes", "circle.xml"))
        elif domain == "L":
            mesh = Mesh(os.path.join(here, "meshes", "l.xml"))
        else:
            raise ValueError, "Unknown domain: {}".format(domain)

        self.boundary_partition = self.boundary_division(domain)

        self.facet_func = FacetFunction("size_t", mesh,
                                        len(self.boundary_partition))

        for k in range(len(self.boundary_partition)):
            self.boundary_partition[k].mark(self.facet_func, k)

        return mesh

    def boundary_division(self, domain):
        """
        Takes as argument a mesh name, and returns
        a partition of it into subdomains.
        """

        class bdy00(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0], 0) and on_boundary
        class bdy01(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0], 1) and on_boundary

        class bdy10(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1], 0) and on_boundary
        class bdy11(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[1], 1) and on_boundary

        class bdy20(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[2], 0) and on_boundary
        class bdy21(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[2], 1) and on_boundary

        class dolfin_interior(SubDomain):
            def inside(self, x, on_boundary):
                return (not (near(x[0], 0) or near(x[0], 1) or
                            near(x[1], 0) or near(x[1], 1))
                        and on_boundary)

        # for domains we don't specify a partition for,
        # use the entire domain
        class entire_boundary(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary
        bdy_partition = {
            "UnitInterval": [bdy00(), bdy01()],
            "UnitSquare": [bdy00(), bdy01(), bdy10(), bdy11()],
            "UnitCube": [bdy00(), bdy01(), bdy10(), bdy11(), bdy20(), bdy21()],
            "Dolfin": [bdy00(), bdy01(), bdy10(), bdy11(),
                       dolfin_interior()]
        }

        if domain in bdy_partition:
            return bdy_partition[domain]
        else:
            return [entire_boundary()]


    def get_bcs(self, V, default="0"):
        """
        Creates Dirichlet boundary conditions using the values of
        bdyK given in self.params, and returns a list with a BC for
        each piece of the boundary. V is the FunctionSpace
        in which the solution lives.
        """
        bcs = []
        for k in range(len(self.boundary_partition)):
            try:
                bc_k_str = self.params["bdy{}".format(k)]
                if bc_k_str.lower() is "none":
                    continue
                bc_expr = self.s2d(bc_k_str)
            except:
                bc_expr = self.s2d(default)
            bcs.append(DirichletBC(V, bc_expr, self.facet_func, k))
        return bcs


    def plot(self):
        # Plot solution
        tmpfile = tempfile.NamedTemporaryFile(dir='/tmp', delete=False,
                  suffix=".png", prefix="fenicsbot_")
        tmpfile_name = tmpfile.name
        tmpfile.close()

        try:
            parameters["plotting_backend"] = "matplotlib"
        except:
            pass
        plot(self.solution)
        plt.savefig(tmpfile_name[:-4], bbox_inches="tight")
        plt.close()

        return tmpfile_name
