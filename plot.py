from types import NoneType

from result import Result
import matplotlib.pyplot as plt

path = 'Results/qt/with_qt'
data = Result.load_results(path)
data.slip_rate_evolution(save = False, path = 'images/Taux/', name='article') # Plot slip rate evolution58479002277847_1.0.pkl") # Load previously saved results in the Res
data.displacement_evolution(save = False, path = 'images/Taux/', name='article') # Plot displacement evolution
data.phase_portrait(save = False, path = 'images/complete/', name='article') # Plot phase portrait
data.theta_evolution(save = False, path = 'images/complete/') # Plot state variable evolution
if data.Tau.__class__ != NoneType:
    data.tau_evolution(save = False, path = 'images/complete/') # Plot shear stress evolution
if data.Sigma_n.__class__ != NoneType:
    data.sigma_evolution(save = False, path = 'images/complete/') # Plot normal stress evolution
if data.P.__class__ != NoneType:
    data.pressure_evolution(save = False, path = 'images/complete/', name='article') # Plot pressure evolution
if data.q.__class__ != NoneType:
    data.flow_rate(save = False, path = 'images/complete/', name='article')
plt.show()

data.magnitude() # Compute moment magnitude
data.cycles() # Compute frequencies of slip events