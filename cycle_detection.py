from result import Result
import numpy as np
import matplotlib.pyplot as plt

path = 'Results/PR_QDYN_RNS/tests.pkl'
data = Result.load_results(path)
taux=0.5
nb_periode_pour_securite=1

signal = data.Phi
time = data.T

period = data.fft_dominant_period()

N=len(time)
max=signal.max()
min=signal.min()
threshhold=taux*(max-min)+min
duree_securite=nb_periode_pour_securite*period
temps_debut=duree_securite
temps_fin=time[-1]-duree_securite

indice_debut=0
indice_fin=N-1

while time[indice_debut]<=temps_debut:
    indice_debut+=1

while time[indice_fin]>=temps_fin:
    indice_fin-=1

indice_debut_injection=np.argmin(np.abs(signal[indice_debut:indice_fin]-threshhold))+indice_debut

plt.plot(time,signal,color='k')
plt.scatter(time[indice_debut_injection],signal[indice_debut_injection],marker='+',s=1000)
plt.grid()
plt.show()