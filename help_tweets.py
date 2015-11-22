DOCUMENTATION_MESSAGE = "You can find a description of the things I can help you with at http://bit.ly/1T4KuNt."    
DOMAIN_MESSAGE = "The domains I know about are UnitInterval, UnitSquare, UnitCube, Dolfin and Circle. Specify them with domain=<domain name>."
BC_MESSAGE = "To specify a boundary condition, use the syntax bdyN=<expression>. I currently only support Dirichlet BCs."
POISSON_MESSAGE = "The Poisson equation describes heat conduction (and other things). When you ask me to solve it, you can specify the forcing term f."
LAPLACE_MESSAGE = "The Laplace equation is a special case of the Poisson equation where f=0."


help_dict = {
    "poisson": POISSON_MESSAGE,
    "laplace": LAPLACE_MESSAGE,
    "domain": DOMAIN_MESSAGE,
    "documentation": DOCUMENTATION_MESSAGE,
    "doc": DOCUMENTATION_MESSAGE,
    "docs": DOCUMENTATION_MESSAGE,
    "boundary": BC_MESSAGE,
    "bc": BC_MESSAGE,
    "bcs": BC_MESSAGE
}

if __name__ == "__main__":
    # Verify that all messages are short enough to be tweeted.
    for subject in help_dict:
        if len(help_dict[subject]) > 140:
            print ("Error: Help message for {} is\n{}\n"
                   "which is too long.").format(subject, help_dict[subject])
