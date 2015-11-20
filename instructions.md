The Twitter FEniCS bot
==============
[Visit the Twitter FEniCs bot](https://twitter.com/fenicsbot/)

With FEniCs bot you can tweet your PDEs problems [@fenicsbot](https://twitter.com/fenicsbot/) and you will receive the solution, obtained with FEniCS, displayed as an answer!

Syntax:
```
@fenicsbot Solve _Problem_ with domain=_Domain_ and f=_ExternalForce_
```

Valid inputs:

`_Problem_`:

1. `Poisson`
2. `Stokes`
3. `LinearElasticity`

`_Domain_`:

1. `UnitInterval`
2. `UnitSquare`
3. `UnitCube`
4. `Dolfin`

`_ExternalForce_`:
For Poisson you need to pass a forcing scalar, such as:

1. Expression: `f=pi*cos(pi*x[0])*cos(pi*x[1])`
2. Constant: `f=0.0`

For Stokes and LinearElasticity you need to pass a forcing vector, such as:

1. Expression: `f=cos(pi*x[0]),sin(pi*x[1])`
2. Constant: `f=0,1`


Example
-------
Poisson equation
```
@fenicsbot Solve Poisson with domain=UnitSquare and f=pi*cos(pi*x[0])*cos(pi*x[1])
```

Stokes equations
```
@fenicsbot Solve Stokes with domain=UnitSquare and f=pi*cos(pi*x[0])*cos(pi*x[1]),10
```


Made by [Karl Erik Holter](https://twitter.com/karl__erik), [Eleonora Piersanti](https://twitter.com/eleonorapiersan) and [Simon Funke](https://twitter.com/SimonFunke) at the BioComp Simula Hackathon 2015 at Finse, Norway.

