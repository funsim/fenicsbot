from solvers import *

# this should be imported from solvers
solvers_by_name = {
    "Poisson": PoissonSolver,
    "Stokes": StokesSolver,
    "LinearElasticity": LinearElasticitySolver
}


def parse(tweet):
    """
    :param tweet: tweet to parse - expected to be in the format
              "@FEniCSbot Solve <SolverName> with <par1>=p1 and <par2>=p2 and ..."
    """
    tweet = excise(tweet)

    try:
        solver_name = tweet.split(" with ")[0]
        specified_params = tweet.split(" with ")[1].split(" and ")
    except:
        # if this happens, the tweet should be just 
        # "@fenicsbot Solve <SolverName>", and we've already excised the preamble
        solver_name = tweet
        specified_params = {}

    solver_name = solver_name.strip()

    
    solver = solvers_by_name[solver_name]

    param_dict = solver.default_parameters()
    
    for p in specified_params:
        parname, parval = p.strip().split("=")
        param_dict[parname] = parval


    solver = solver(param_dict)
    
    return solver

def excise(s):
    # print s
    username = "@fenicsbot solve "
    s_lower = s.lower()
    # print s_lower
    i = s_lower.find(username)
    # print s, i
    # assert i > -1
    s = s[:i] + s[i+len(username):]
    # print s
    return s


