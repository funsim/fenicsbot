""" A simple test script to simulate an handling request """

from parser import parse

tweet = "@fenicsbot solve Poisson in domain UnitSquare with f=1/(x+0.5)"
solver = parse(tweet)
solver.solve()
print solver.plot()
