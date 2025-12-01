from result1D import Result1D
import matplotlib.pyplot as plt

path = 'model2/Results1D/test/test_simulation.pkl'

data = Result1D.load_results(path)
data.theta_evolution() # Plot state variable evolution

plt.show()