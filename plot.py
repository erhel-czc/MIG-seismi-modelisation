from result import Result
import matplotlib.pyplot as plt

data = Result.load_results("0.6_1e-11_0.4.pkl") # Load previously saved results in the Results folder
data.slip_rate_evolution() # Plot slip rate evolution
data.phase_portrait() # Plot phase portrait
plt.show()