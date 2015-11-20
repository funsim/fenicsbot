from demos import *



demos_by_name = {
    "Poisson": PoissonSolver,
    "Stokes": StokesSolver
}

def parser(s):
    """
    :param s: tweet to parse - expected to be in the format
              "@FEniCSbot Solve <DemoName> with <par1>=p1; <par2>=p2;..."
    """
    s = excise(s)

    # demo_name, specified_params = s.split("; ")
    try:
        demo_name = s.split(" with ")[0]
        specified_params = s.split(" with ")[1].split(" and ")
    except:
        # if this happens, the tweet should be just 
        # @fenicsbot Solve <DemoName>, and we've already excised
        demo_name = s
        specified_params = {}

    demo_name = demo_name.strip()



    
    demo = demos_by_name[demo_name]

    param_dict = demo.default_parameters()
    
    # print specified_params
    for p in specified_params:
        # print p
        parname, parval = p.strip().split("=")
        param_dict[parname] = parval
        # print parname, param_dict[parname]

    solver = demo(param_dict)
    solver.solve()
    return solver.plot()


     

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

if __name__=="__main__":
    tweet = "@FEniCSbot Poisson2D"
    # tweet = "@FEniCSbot Poisson; D=1; f=1"
    # img_fn = parser(tweet)
    # print img_fn
    
    tweets = [
        "@fenicsbot Solve Poisson with f=1",
        "@fenicsbot Solve Poisson",
        "@fenicsbot Solve Poisson with f=1 and D=2",
        "@fenicsbot Solve Stokes"
    ]


    for t in tweets:
        parser(t)
    
    # for s, s_e in [
    #         ("@fenicsbot Solve text", "text"),
    #         ("text @fenicsbot Solve text", "text text")]:
    #     # print s_e
    #     # print s
    #     # print excise(s)
    #     assert s_e == excise(s)
