import matplotlib.pyplot as plt
import seaborn
import numpy as np

data = np.genfromtxt("C_HB.dat")
x = data[:, 0]
y = data[:, 1]
yerr = data[:, 2]
plt.errorbar(x, y, yerr=yerr, fmt="x")
plt.scatter(x, y, marker="x", label="heatbath")

data = np.genfromtxt("C_MP.dat")
x = data[:, 0]
y = data[:, 1]
yerr = data[:, 2]
plt.errorbar(x, y, yerr=yerr, fmt="x")
plt.scatter(x, y, marker="x", label="Metropolis")
plt.xlim(
    0,
)
plt.xlabel("T")
plt.ylabel("C")

plt.legend()
plt.show()
