from result import Result
import numpy as np
import matplotlib.pyplot as plt

def find_cycle_start(path_without_pressure, taux):
    data = Result.load_results(path_without_pressure)
    nb_periode_pour_securite_debut=1
    nb_periode_pour_securite_fin=3
    charge=False

    signal = np.array(data.Phi)
    time = data.T
    period = data.fft_cycle_period()

    N=len(time)
    max=signal.max()
    min=signal.min()
    threshhold=taux*(max-min)+min
    temps_debut=nb_periode_pour_securite_debut*period
    temps_fin=time[-1]-nb_periode_pour_securite_fin*period

    indice_debut=0
    indice_fin=N-1

    while time[indice_debut]<=temps_debut:
        indice_debut+=1

    while time[indice_fin]>=temps_fin:
        indice_fin-=1

    decalage_debut=indice_debut 
    signal_tronque=signal[indice_debut:indice_fin].copy()
    temps_tronque=time[indice_debut:indice_fin].copy()

    if charge == True:
        mask= np.diff(signal_tronque)>=0
        signal_tronque[1:][np.logical_not(mask)]=np.inf
        indice_debut_injection=(np.abs((signal_tronque[+1:]-threshhold))).argmin()
        indice_debut_injection+=decalage_debut+1

    if charge == False:
        mask= np.diff(signal_tronque)<=0
        signal_tronque[1:][np.logical_not(mask)]=np.inf
        indice_debut_injection=(np.abs((signal_tronque[+1:]-threshhold))).argmin()
        indice_debut_injection+=decalage_debut+1

    plt.figure(2)
    plt.plot(time,signal,color='k')
    plt.scatter(time[indice_debut_injection],signal[indice_debut_injection],marker='+',s=100)
    plt.grid()
   #plt.show()

    return (time[indice_debut_injection]+time[indice_debut_injection-1])/2