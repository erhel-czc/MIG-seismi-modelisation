from result import Result
from matplotlib import pyplot as plt


for i in range(10):
    path=f"Results/taux_0_to_1/Taux=0.{i}0"
    data = Result.load_results(path)
    cycles=Result.cycles(data)
    plt.plot(cycles,'+')

""""
path=f"Results/taux_0_to_1/Taux=0.10"
data = Result.load_results(path)
cycles1=Result.cycles(data)

path=f"Results/taux_0_to_1/Taux=0.20"
data = Result.load_results(path)
cycles2=Result.cycles(data)

path=f"Results/taux_0_to_1/Taux=0.30"
data = Result.load_results(path)
cycles3=Result.cycles(data)

path=f"Results/taux_0_to_1/Taux=0.40"
data = Result.load_results(path)
cycles4=Result.cycles(data)

path=f"Results/taux_0_to_1/Taux=0.50"
data = Result.load_results(path)
cycles5=Result.cycles(data)

path=f"Results/taux_0_to_1/Taux=0.20"
data = Result.load_results(path)
cycles6=Result.cycles(data)

path=f"Results/taux_0_to_1/Taux=0.70"
data = Result.load_results(path)
cycles7=Result.cycles(data)

path=f"Results/taux_0_to_1/Taux=0.80"
data = Result.load_results(path)
cycles8=Result.cycles(data)

path=f"Results/taux_0_to_1/Taux=0.90"
data = Result.load_results(path)
cycles9=Result.cycles(data)

path=f"Results/taux_0_to_1/Taux=1.00"
data = Result.load_results(path)
cycles10=Result.cycles(data)

plt.plot(cycles1,'+',label=f"taux = 0.1")
plt.plot(cycles2,'+',label=f"taux = 0.2")
plt.plot(cycles3,'+',label=f"taux = 0.3")
plt.plot(cycles4,'+',label=f"taux = 0.4")
plt.plot(cycles5,'+',label=f"taux = 0.5")
plt.plot(cycles6,'+',label=f"taux = 0.6")
plt.plot(cycles7,'+',label=f"taux = 0.7")
plt.plot(cycles8,'+',label=f"taux = 0.8")
plt.plot(cycles9,'+',label=f"taux = 0.9")
plt.plot(cycles10,'+',label=f"taux = 1.0")
"""

plt.show()