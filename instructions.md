The Twitter FEniCS bot
==============
[Visit the Twitter FEniCs bot](https://twitter.com/fenicsbot/)

With FEniCs bot you can tweet your PDEs problems [@fenicsbot](https://twitter.com/fenicsbot/) and you will receive the solution, obtained with FEniCS, displayed as an answer!

Tweet Syntax
------------
```
@fenicsbot Solve _Problem_ with par1=_par1val_ and par2=_par2val_ and ...
```

The valid options are problem specific.


Valid inputs:

`_Problem_`:

1. `Poisson`, options: `f`: external force, `domain`: domain, `bdyK`: Dirichlet BC on piece K of the boundary
2. `Stokes`, options: `f`: external force, `domain`: domain
3. `LinearElasticity`: `f`: external force, `domain`: domain, `E`: Youngs modulus, `nu`: Poissons ratio
4. `Burgers`: `f`: external force, `domain`: domain, `ic`: initial condition, `dt`: timestep, `T`: final time, `nu`: viscosity

`_Domain_`:

1. `UnitInterval`
2. `UnitSquare`
3. `UnitCube`
4. `Dolfin`
4. `Circle`
4. `L`

`_ExternalForce_`:
For Poisson you need to pass a forcing scalar, such as:

1. Expression: `f=pi*cos(pi*x[0])*cos(pi*x[1])`
2. Constant: `f=0.0`

For Stokes and LinearElasticity you need to pass a forcing vector, such as:

1. Expression: `f=cos(pi*x[0]),sin(pi*x[1])`
2. Constant: `f=0,1`

`_BCs_`:
For solvers including `bdy0`, `bdy1`, ... as options, you can impose Dirichlet boundary conditions on reasonable pieces of the boundary. In order to impose the condition u=g on piece number 2 of the boundary, pass the option `bdy2=g`. All BCs not given default to 0.

Example
-------
Poisson equation
```
@fenicsbot Solve Poisson with domain=UnitSquare and f=pi*cos(pi*x[0])*cos(pi*x[1])
```
```
@fenicsbot Solve Poisson with domain=UnitSquare and f=0 and bdy0=1 and bdy1=1
```

Stokes equations
```
@fenicsbot Solve Stokes with domain=UnitSquare and f=pi*cos(pi*x[0])*cos(pi*x[1]),10
```

LinearElasticity
```
@fenicsbot Solve LinearElasticity with domain=Dolfin and f=sin(pi*x[0]),cos(pi*x[1])
```

Made by [Karl Erik Holter](https://twitter.com/karl__erik), [Eleonora Piersanti](https://twitter.com/eleonorapiersan) and [Simon Funke](https://twitter.com/SimonFunke) at the BioComp Simula Hackathon 2015 at Finse, Norway.

