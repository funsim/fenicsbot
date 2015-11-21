""" A simple test script to simulate an handling request """

import sys
from parser import parse

if __name__ == "__main__":
    tweet = sys.argv[1]
    #tweet = "@fenicsbot solve Poisson in domain UnitSquare with f=1/(x+0.5)"
    solver = parse(tweet)
    solver.solve()
    print solver.plot()
