import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import os
import tempfile
from dolfin import *
from mshr import UnitSphereMesh
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

    def get_mesh(self, return_bdys=False):
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

        bdy_subdiv = self.boundary_division(domain)
        bdy_enum = FacetFunction("size_t", mesh, len(bdy_subdiv))

        for k in range(len(bdy_subdiv)):
            bdy_subdiv[k].mark(bdy_enum, k)

        if return_bdys:
            return mesh, bdy_subdiv, bdy_enum
        else:
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
            return []




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
