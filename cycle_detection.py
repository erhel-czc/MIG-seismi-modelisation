from result import Result
import numpy as np
import matplotlib.pyplot as plt

path = 'Results/PR_QDYN_RNS/tests.pkl'
data = Result.load_results(path)

data.slip_rate_evolution()

data.fft_cycle_period(output=True)