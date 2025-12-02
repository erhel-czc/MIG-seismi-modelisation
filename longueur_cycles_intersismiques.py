from result import Result
from matplotlib import pyplot as plt

plt.figure("Longueur des cycles intersismiques")

for i in range(10):
    path=f"Results/taux_0_to_1/Taux=0.{i}0"
    data = Result.load_results(path)
    cycles=Result.cycles(data)
    plt.plot(cycles,'+', label=f"taux = 0.{i}0", markersize=10)

plt.title("Longueur des cycles intersismiques en fonction du numéro de cycle")
plt.xlabel("Numéro de cycle")
plt.ylabel("Longueur des cycles intersismiques")
plt.legend()
plt.grid()

plt.savefig("images/longueur_cycles_intersismiques.pdf")

plt.show()