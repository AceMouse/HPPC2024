import numpy as np
import matplotlib.pyplot as plt

# loading data
data = np.loadtxt("weak_scale.dat", delimiter=' ', skiprows=1)
n = data[:,0]
total = data[:,1]
fft = data[:,2]
non_fft = data[:,3]

# calculating speedup
fft_su = fft[0]/fft
non_fft_su = non_fft[0]/non_fft
total_su = total[0] / total


# plotting weak scaling
fig,ax = plt.subplots()
ax.plot(n, total_su, label = "Total")
ax.plot(n, fft_su, label = "FFT")
ax.plot(n, non_fft_su, label = "Non FFT")
ax.set_title("Efficiency (speedup with increasing problem size)")
ax.set_xlabel("Number of threads")
ax.set_xticks(n[1:])
ax.set_ylabel("Efficiency")
ax.legend()
plt.tight_layout()
plt.savefig("./img/weak_scaling.png")


# loading data
data = np.loadtxt("strong_scale.dat", delimiter=' ', skiprows=1)
s_total = data[:,1]
s_fft = data[:,2]
s_non_fft = data[:,3]

# calculating speedup
s_total_su = s_total[0] / s_total

# plotting strong scaling
fig,ax = plt.subplots()
ax.plot(n, s_total_su)
ax.set_title("Speedup for strong scaling")
ax.set_xlabel("Number of threads")
ax.set_xticks(n[1:])
ax.set_ylabel("Speedup")
plt.tight_layout()
plt.savefig("./img/strong_scaling.png")


# calc parallel frac
p_frac = (1/s_total_su[-1] - 1)/(1/n[-1] -1)
print(f"Speedup with {n[-1]} cores is {total_su[-1]:.3f}")
print(f"Parallel fraction {p_frac:.4f}, serial {1-p_frac:.4f}")

# plotting weak scaling efficiency
ef_total = s_total_su/n
fig,ax = plt.subplots()
ax.plot(n, ef_total)
ax.set_title("Efficiency for strong scaling")
ax.set_xlabel("Number of threads")
ax.set_ylim(0, 1.05)
ax.set_xticks(n[1:])
ax.set_ylabel("Efficiency")
plt.tight_layout()
plt.savefig("./img/effi_strong_scaling.png")

