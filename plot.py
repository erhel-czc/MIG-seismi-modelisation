from result import Result
import matplotlib.pyplot as plt

path = '/Users/rico/Desktop/Mines/MIG/Semaine 2/MIG-seismi-modelisation/Results/PR_QDYN_RNS_modele_oriente/0.5_1.0339460065412583e-08_0.01258479002277847_1.0.pkl'
data = Result.load_results(path)
data.slip_rate_evolution() # Plot slip rate evolution58479002277847_1.0.pkl") # Load previously saved results in the Res
data.phase_portrait() # Plot phase portrait
data.tau_evolution() # Plot shear stress evolution
data.sigma_evolution() # Plot normal stress evolution
plt.show()