from result import Result
import matplotlib.pyplot as plt
import numpy as np


def plot_Mw_vs_lenght_fault(num_save):
    paths = [f'Results/PR_QDYN_RNS_modele_oriente/Lenght_Fault_varies/Lenght Fault {i}' for i in range(num_save)]
    data = [Result.load_results(path) for path in paths]
    lenght_faults = np.array([d.pd.lenght_fault for d in data])
    Mw = np.array([d.magnitude(print_result=False) for d in data])

    plt.figure('Mw vs Lenght Fault')

    plt.scatter(lenght_faults, Mw, color='black', marker='+')
    plt.plot(lenght_faults, Mw,'-k', alpha=0.3)
    plt.xlabel('Lenght Fault (m)')
    plt.ylabel('Moment Magnitude (Mw)')
    plt.title('Moment Magnitude vs Lenght Fault')
    plt.grid()
    
    plt.show()

plot_Mw_vs_lenght_fault(100)