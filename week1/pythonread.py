import numpy as np
import matplotlib.pyplot as plt
x, y, z = np.loadtxt('sir_output.txt', skiprows=1, unpack=True)

step_size = 0.5
days = np.arange(0, len(x)*step_size, step_size)

plt.plot(days, x, label = "S")
plt.plot(days, y, label = "I")
plt.plot(days, z, label = "R")
plt.title("SIR model")
plt.ylabel("Individuals")
plt.xlabel("Days (since first infection)")
plt.legend()
plt.savefig("sir_plot.png")