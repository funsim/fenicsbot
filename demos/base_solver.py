import os
import tempfile
from dolfin import *

class BaseSolver(object):

    @staticmethod
    def default_parameters():
        raise NotImplementedError

    def solve(self):
        raise NotImplementedError

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
        else:
            raise ValueError, "Unknown domain: {}".format(domain)

        return mesh

    def plot(self):
	# Plot solution
        tmpfile = tempfile.NamedTemporaryFile(dir='/tmp', delete=False,
                  suffix=".png", prefix="fenicsbot_")
        tmpfile_name = tmpfile.name
        tmpfile.close()

	plot(self.solution).write_png(tmpfile_name[:-4])
        return tmpfile_name
