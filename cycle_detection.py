# a program containing cycle detection functions
"""
CONSIDÉRER QU'À PARTIR DE LA MOITIÉ DU TEMPS
"""

from result import Result
import numpy as np
import matplotlib.pyplot as plt

path = 'Results/PR_QDYN_RNS/tests.pkl'
data = Result.load_results(path)
data.slip_rate_evolution()


signal = data.V
time = data.T

sampling_rate = 1 / np.mean(np.diff(time))
fft = np.fft.rfft(signal)
freqs = np.fft.rfftfreq(len(signal), d=1/sampling_rate)
periode = 1 / freqs[np.argmax(np.abs(fft[1:])) + 1]
print(periode)

plt.figure('FFT of v(t)')
plt.plot(freqs[1:], np.abs(fft[1:]))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('FFT of v(t)')
plt.grid()

plt.show()