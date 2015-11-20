The FEniCS bot
==============
[FEniCs bot](https://twitter.com/fenicsbot/)

With FEniCs bot you can tweet your PDEs problems [@fenicsbot](https://twitter.com/fenicsbot/) and you will receive the solution, obtained with FEniCS, displayed as an answer!

Syntax:
```
@fenicsbot Solve _Problem_ with domain=_Domain_ and f=_ExternalForce_
```

_Problem_:

1. Poisson 
2. Stokes 
3. Linear Elasticity 

_Domain_:

1. UnitInterval
2. UnitSquare
3. UnitCube
4. Dolfin

_ExternalForce_:

1. `Expression: f = pi*cos(pi*x[0])*cos(pi*x[1])`
2. `Constant: f = 0.0`

Example:
@fenicsbot Solve Poisson with domain=UnitSquare and `f=pi*cos(pi*x[0])*cos(pi*x[1])`


