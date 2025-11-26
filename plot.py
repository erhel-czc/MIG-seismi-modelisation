from PR_QDYN_RNS import Result, ParamMec, NdParamMec, ParamComp

data = Result.load_results("testPressionv2.pkl") #Load previously saved results in the Results folder
data.slip_rate_evolution()  #Plot slip rate evolution
data.phase_portrait()      #Plot phase portrait