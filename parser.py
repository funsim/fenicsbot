from demos import *


demos_by_name = {
    "Poisson": PoissonSolver,
    "Stokes": StokesSolver
}


def parse(tweet):
    """
    :param tweet: tweet to parse - expected to be in the format
              "@FEniCSbot Solve <DemoName> with <par1>=p1 and <par2>=p2 and ..."
    """
    tweet = excise(tweet)

    try:
        demo_name = tweet.split(" with ")[0]
        specified_params = tweet.split(" with ")[1].split(" and ")
    except:
        # if this happens, the tweet should be just 
        # "@fenicsbot Solve <DemoName>", and we've already excised the preamble
        demo_name = tweet
        specified_params = {}

    demo_name = demo_name.strip()

    
    demo = demos_by_name[demo_name]

    param_dict = demo.default_parameters()
    
    for p in specified_params:
        parname, parval = p.strip().split("=")
        param_dict[parname] = parval


    solver = demo(param_dict)
    
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


