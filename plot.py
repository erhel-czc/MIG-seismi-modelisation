from result import Result
import matplotlib.pyplot as plt

path = 'Results/PR_QDYN_RNS_P/pressure1'
data = Result.load_results(path)
data.slip_rate_evolution(save = True, path = 'images/pressure/', name='pressure1') # Plot slip rate evolution58479002277847_1.0.pkl") # Load previously saved results in the Res
data.phase_portrait(save = True, path = 'images/pressure/', name='pressure1') # Plot phase portrait
if data.Tau != None:
    data.tau_evolution() # Plot shear stress evolution
if data.Sigma_n != None:
    data.sigma_evolution() # Plot normal stress evolution
plt.show()

data.magnitude() # Compute moment magnitude