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



    def get_mesh(self):
        domain = self.params["domain"]

        if domain == "UnitInterval":
                mesh = UnitIntervalMesh(20)
        elif domain == "UnitSquare":
                   mesh = UnitSquareMesh(20, 20)
        elif domain == "UnitCube":
                mesh = UnitCubeMesh(5, 5, 5)
        elif domain == "Dolfin":
                here = os.path.dirname(__file__)
                mesh = Mesh(os.path.join(here, "dolfin.xml.gz"))
        elif domain == "Sphere":
                mesh = UnitSphereMesh(10)
                plot(mesh, interactive=True)
        else:
            raise ValueError, "Unknown domain: {}".format(domain)

        return mesh

    def plot(self):
        # Plot solution
        tmpfile = tempfile.NamedTemporaryFile(dir='/tmp', delete=False,
                  suffix=".png", prefix="fenicsbot_")
        tmpfile_name = tmpfile.name
        tmpfile.close()

        parameters["plotting_backend"] = "matplotlib"
        plot(self.solution)
        plt.savefig(tmpfile_name[:-4], bbox_inches="tight")

        return tmpfile_name
