# a program containing cycle detection functions
"""
CONSIDÉRER QU'À PARTIR DE LA MOITIÉ DU TEMPS
"""

from result import Result
import numpy as np
import matplotlib.pyplot as plt

path = 'Results/PR_QDYN_RNS/tests.pkl'
data = Result.load_results(path)
# data.slip_rate_evolution()
taux=0.5
nb_periode_pour_securite=1

signal = data.Phi
time = data.T

sampling_rate = 1 / np.mean(np.diff(time))
fft = np.fft.rfft(signal)
freqs = np.fft.rfftfreq(len(signal), d=1/sampling_rate)
periode = 1 / freqs[np.argmax(np.abs(fft[1:])) + 1]

"""""
print(periode)


plt.figure('FFT of v(t)')
plt.plot(freqs[1:], np.abs(fft[1:]))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('FFT of v(t)')
plt.grid()

plt.show()
"""""


N=len(time)
max=signal.max()
min=signal.min()
threshhold=taux*(max-min)+min
duree_securite=nb_periode_pour_securite*periode
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
=======

data.slip_rate_evolution()

data.fft_cycle_period(output=True)
>>>>>>> bc86bed1d641a70c8a76c9cd9dee7a7c9560ea7e
