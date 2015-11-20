""" A simple test script to simulate an handling request """

from parser import parse

tweet = "@fenicsbot Solve Poisson with f=sin(x[0])*sin(x[1]) and domain=Dolfin"
solver = parse(tweet)
solver.solve()
solver.plot()
