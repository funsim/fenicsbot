import re
from solvers import *

# this should be imported from solvers
solvers_by_name = {
    "Poisson": PoissonSolver,
    "Stokes": StokesSolver,
    "LinearElasticity": LinearElasticitySolver,
    "Burgers": BurgersSolver,
}


def parse(tweet):
    """
    :param tweet: tweet to parse - expected to be in the format
              "@FEniCSbot Solve <SolverName> with <par1>=p1 and <par2>=p2 and ..."
    """
    tweet = excise(tweet)

    # do some replacement to handle some common cases
    tweet = tweet.replace("'", "")
    tweet = tweet.replace('"', "")
    tweet = tweet.replace(" and ", " with ")
    tweet = tweet.replace(" in ", " with ")
    tweet = tweet.replace(" on ", " with ")
    # Remove whitespaces before and after = signs
    tweet = re.sub('\s*=\s*', '=', tweet)

    tweet = tweet.replace("domain ", "domain=")

    try:
        tweetl = tweet.split(" with ")
        solver_name, specified_params = tweetl[0], tweetl[1:]
    except:
        # if this happens, the tweet should be just
        # "@fenicsbot Solve <SolverName>", and we've already excised the preamble
        solver_name = tweet
        specified_params = {}

    solver_name = solver_name.strip()

    try:
        solver = solvers_by_name[solver_name]
    except:
        raise ValueError, "{} not implemented (yet!)".format(solver_name)

    param_dict = solver.default_parameters()

    for p in specified_params:
        # strip() is to not have trouble with spaces - should allow something
        # like "f = expr" as well as "f= expr", "f =expr" or "f=expr"
        parname, parval = p.strip().split("=")
        parname = parname.strip()
        parval = parval.strip()
        param_dict[parname] = parval


    solver = solver(param_dict)

    return solver

def excise(s):
    username = "@fenicsbot solve "
    i = s.lower().find(username)
    s = s[:i] + s[i+len(username):]
    return s


if __name__ == "__main__":
    error_tweet = "@fenicsbot Solve Poisson with f=sin(x[0])*sin(x[1]) and domain=Dolfin"
    solver = parse(error_tweet)
    solver.solve()
    print solver.plot()
