import pytest
import numpy as np

from fenicsbot.demos import *

@pytest.mark.parametrize("mesh", [
    "UnitInterval",
    "UnitSquare",
    "UnitCube",
    "Dolfin"
])
def test_poisson_zero_solution(mesh):
    """
    Tests that a Laplace equation with identically zero Dirichlet 
    boundary conditions has solution u=0.
    """
    
    solver = PoissonSolver({"mesh": mesh, "f": 0})
    solver.solve()
    u = solver.solution.vector().array()
    
    assert np.allclose(np.zeros(u.shape), u)


    
