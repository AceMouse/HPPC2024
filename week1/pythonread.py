import numpy as np
import matplotlib.pyplot as plt
x, y, z = np.loadtxt('sir_output.txt', skiprows=1, unpack=True)
plt.plot(x, label = "S")
plt.plot(y, label = "I")
plt.plot(z, label = "R")
plt.title("SIR model")
plt.ylabel("Individuals")
plt.xlabel("Days (since first infection)")
plt.legend()
plt.savefig("sir_plot.png")