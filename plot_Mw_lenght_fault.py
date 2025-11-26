from result import Result
import matplotlib.pyplot as plt
import numpy as np

n = 100

def load_all_results(num_save):
    paths = [f'Results/PR_QDYN_RNS_modele_oriente/Lenght_Fault_varies/Lenght Fault {i}' for i in range(num_save)]
    data = [Result.load_results(path) for path in paths]
    return data

def plot_Mw_vs_lenght_fault(data, save=False):
    lenght_faults = np.array([d.pd.lenght_fault for d in data])
    Mw = np.array([d.magnitude(print_result=False) for d in data])

    plt.figure('Mw vs Lenght Fault')

    plt.scatter(lenght_faults, Mw, color='black', marker='+')
    plt.plot(lenght_faults, Mw,'-k', alpha=0.3)
    plt.xlabel('Lenght Fault (m)')
    plt.ylabel('Moment Magnitude (Mw)')
    plt.title('Moment Magnitude vs Lenght Fault')
    plt.grid()
    
    if save : plt.savefig('images/Mw_vs_lenght_fault.pdf')

def plot_stress_drop_vs_lenght_fault(data, save=False):
    # stress_drop*L**3=M0
    lenght_faults = np.array([d.pd.lenght_fault for d in data])
    Mw = np.array([d.magnitude(print_result=False) for d in data])
    
    # Filter out None values
    valid_idx = np.array([i for i, m in enumerate(Mw) if m is not None])
    lenght_faults = lenght_faults[valid_idx]
    Mw = Mw[valid_idx]
    
    M0 = 10**(1.5*(Mw + 6))
    stress_drops = M0/(lenght_faults**3)

    plt.figure('Stress Drop vs Lenght Fault')

    plt.scatter(lenght_faults, stress_drops, color='black', marker='+')
    plt.plot(lenght_faults, stress_drops,'-k', alpha=0.3)
    plt.xlabel('Lenght Fault (m)')
    plt.ylabel('Stress Drop (Pa)')
    plt.title('Stress Drop vs Lenght Fault')
    plt.grid()
    
    if save : plt.savefig('images/stress_drop_vs_lenght_fault.pdf')

data = load_all_results(n)
plot_Mw_vs_lenght_fault(data)
plot_stress_drop_vs_lenght_fault(data)
#plot_stress_drop_vs_shear(data)
plt.show()