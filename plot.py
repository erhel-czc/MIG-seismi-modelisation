from result import Result
import matplotlib.pyplot as plt

path = 'Results/PR_QDYN_RNS_P/pressure1'
data = Result.load_results(path)
data.slip_rate_evolution() # Plot slip rate evolution58479002277847_1.0.pkl") # Load previously saved results in the Res
data.phase_portrait() # Plot phase portrait
#data.tau_evolution() # Plot shear stress evolution
#data.sigma_evolution() # Plot normal stress evolution
plt.show()