import tempfile
from dolfin import plot

class BaseSolver(object):

    @staticmethod
    def default_parameters():
        raise NotImplementedError

    def solve(self):
        raise NotImplementedError

    def plot(self):
	# Plot solution
        tmpfile = tempfile.NamedTemporaryFile(dir='/tmp', delete=False,
                  suffix=".png", prefix="fenicsbot_")
        tmpfile_name = tmpfile.name
        tmpfile.close()

	plot(self.solution).write_png(tmpfile_name[:-4])
        return tmpfile_name
