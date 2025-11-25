from result import Result
import numpy as np

path = '/Users/rico/Desktop/Mines/MIG/Semaine 2/MIG-seismi-modelisation/Results/PR_QDYN_RNS_modele_oriente/0.5_1.0339460065412583e-08_0.01258479002277847_1.0.pkl'
data = Result.load_results(path) # Load previously saved results in the Results folder

threshold_velocity = 1e-3  # Define a threshold velocity to identify slip events
