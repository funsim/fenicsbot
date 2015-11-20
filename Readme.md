The FEniCS bot
==============

The [@fenicsbot](https://twitter.com/fenicsbot/) listens for tweets with a PDE
problem, solves the PDE and sends the solution as a twitter reply.

Instructions
------------

[See here](instructions.md)


Contribute
----------

Feel free to send you pull requests to make the FEniCS bot even better.

It should be easy to add functionality to the FEniCS bot. For example:

*How to add a new solver*

A new solver can be implemented in two steps:

1. Create a new FEniCS solver in the `solvers` directory (see existing solvers as a template)
2. Open parser.py and add your solver to the `solvers_by_name` dictionary of supported solvers.


*How to test an incoming tweet*

Use `test_tweet.py` to test what FEniCS bot would do for a given twitter message.
