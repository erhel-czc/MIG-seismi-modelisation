from types import NoneType

from result import Result
import matplotlib.pyplot as plt

path = 'Results/PR_QDYN_RNS_P/article_2.50e+06_1.70e+00_1.00e-11_4.00e-01.pkl'
data = Result.load_results(path)
data.slip_rate_evolution(save = True, path = 'images/pressure/', name='article') # Plot slip rate evolution58479002277847_1.0.pkl") # Load previously saved results in the Res
data.phase_portrait(save = True, path = 'images/pressure/', name='article') # Plot phase portrait
if data.Tau.__class__ != NoneType:
    data.tau_evolution() # Plot shear stress evolution
if data.Sigma_n.__class__ != NoneType:
    data.sigma_evolution() # Plot normal stress evolution
if data.P.__class__ != NoneType:
    data.pressure_evolution(save = True, path = 'images/pressure/', name='article') # Plot pressure evolution
plt.show()

data.magnitude() # Compute moment magnitude