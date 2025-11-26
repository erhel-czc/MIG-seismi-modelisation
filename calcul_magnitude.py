from result import Result
import numpy as np

path = '/Users/rico/Desktop/Mines/MIG/Semaine 2/MIG-seismi-modelisation/Results/PR_QDYN_RNS_modele_oriente/0.5_1.0339460065412583e-08_0.01258479002277847_1.0.pkl'
data = Result.load_results(path) # Load previously saved results in the Results folder

threshold_velocity = 1/data.pd.V_p  # Define a threshold velocity to identify slip events

T=data.T
V=data.V
deplacement=0

if max(V) > threshold_velocity:
    while len(T)>0 and V[0] < threshold_velocity:
        T = T[1:]
        V = V[1:]

    while len(T)>0 and V[0] > threshold_velocity:
        deplacement += V[0] * (T[1] - T[0])
        T = T[1:]
        V = V[1:]
    
    deplacement *= data.pd.dc
    print(f"Slip event detected. Total slip displacement: {deplacement}.")   

    def calcul_magnitude(deplacement, pd):
        """
        Calculate the moment magnitude based on slip displacement and physical parameters.

        Parameters:
        deplacement (float): Slip displacement in meters.
        pd (PhysicalData): An instance of the PhysicalData class containing physical parameters.

        Returns:
        float: Calculated moment magnitude.
        """
        mu = pd.shear
        L = pd.lenght_fault

        # Calculate seismic moment (M0)
        M0 = mu * L**2 * deplacement  # in Newton-meters

        # Calculate moment magnitude (Mw)
        Mw = (2/3) * np.log10(M0) - 6.0

        return Mw

    print(f"Calculated moment magnitude: {calcul_magnitude(deplacement, data.pd)}")
else:
    print("No slip event detected.")