from result1D import Result1D
from oneD_CFAULT_RNS import NdParamMec, ParamComp
import matplotlib.pyplot as plt

path = 'model2/Results1D/test/test_simulation_1.pkl'

data = Result1D.load_results(path)
data.theta_evolution(save=True, path='model2/images') # Plot state variable evolution
data.slip_rate_evolution(save=True, path='model2/images') # Plot slip rate evolution
data.slip_rate_evolution_3D(save=True, path='model2/images')
data.slip_rate_evolution_map(save=True, path='model2/images')

plt.show()