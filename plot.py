from PR_QDYN_RNS import Result, ParamMec, NdParamMec, ParamComp
import matplotlib.pyplot as plt

data = Result.load_results("0.6_1e-11_0.41.pkl") # Load previously saved results in the Results folder
data.slip_rate_evolution() # Plot slip rate evolution
data.phase_portrait() # Plot phase portrait
plt.show()