from demos.poisson import *



demos_by_name = {
    "Poisson": PoissonSolver
}

def parser(s):
    """
        
    :param s: tweet to parse 

    """
    # :returns: a tuple (demo, pardict) with the demo class and a complete
    # dictionary of parameters (for non-specified parameters, use 
    # default values from demo class)
   
    # maybe there should be an assert @FEniCSbot found here?
    s = excise(s)
    # print s

    # demo_name, specified_params = s.split("; ")
    demo_name = s.split("; ")[0]
    demo_name = demo_name.strip()
    specified_params = s.split("; ")[1:]
    demo = demos_by_name[demo_name]
    
    param_dict = demo.default_parameters()
    
    # print specified_params
    for p in specified_params:
        # print p
        parname, parval = p.split("=")
        param_dict[parname] = parval
        # print parname, param_dict[parname]
    
    solver = demo(param_dict)
    solver.solve()
    return solver.plot()


     

def excise(s):
    username = "@fenicsbot"
    s_lower = s.lower()
    
    i = s_lower.find(username)
    # assert i > -1
    s = s[:i] + s[i+len(username):]
    return s

if __name__=="__main__":
    tweet = "@FEniCSbot Poisson2D"
    # tweet = "@FEniCSbot Poisson; D=1; f=1"
    # img_fn = parser(tweet)
    # print img_fn
    
    for s, s_e in [
            ("@fenicsbot text", " text"),
            ("text @fenicsbot text", "text  text")]:
        assert s_e == excise(s)
