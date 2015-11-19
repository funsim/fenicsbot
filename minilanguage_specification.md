@FEniCSbot <DemoType>; <DemoPar=expr>; BCs: <BC1> <BC2> ...

<BCi> is specified as {(bdy_part) Dirichlet/Neumann: <expr>}



if a BCi is not specified, it is assumed to be {<BCi> N: 0}



for example, to specify the 2D problem Laplace(u) = 0, with u=1 on x=0, u=0 on x=1 and du/dn = 0 on y=0, 1:


@FEniCSbot Poisson2D; f=0; BCs: {(x=0) D: 1} {(x=1) D: 0} {(y=0) N: 0} {(y=1) N: 0}


(the last two are optional.) 



Alternatively, the 3D problem Laplace(u) = -xyz with BCs u=x*y on z=0, z=1, u=2y on x=1 and on all other sides du/dn = 0 would be specified as

@FEniCSbot Poisson3D; f=-x*y*z; BCs: {(z=0) D: x*y} {(z=1) D: x*y} {(x=1) N: 2*y}


Ideas:
--------
- syntax like <cont> at the end of a tweet if one runs out of space - i.e. if the tweet

@FEniCSbot text1 text2 text3

is over 140 characters, one could instead first tweet

@FEniCSbot text1 <cont>

then

@FEniCSbot text2 <cont>

then

@FEniCSbot text3

and FEniCSbot would DWYM.
