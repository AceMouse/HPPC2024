import numpy as np
import matplotlib.pyplot as plt

# loading data
data = np.loadtxt("strong_v256.dat", delimiter=' ', skiprows=1)
n = data[:,0]
par = data[:,1]
apar = data[:,2]

# calculating speedup
par_su = par[0]/par
apar_su = apar[0]/apar

# plotting weak scaling
fig,ax = plt.subplots()
ax.plot(n, par_su, label = "Parallelized")
ax.plot(n, apar_su, label = "Par + Async")
ax.set_title("Strong scaling")
ax.set_xlabel("Number of gangs")
ax.set_xticks(n[2:])
ax.set_ylabel("Speedup")
ax.legend()
ax.set_ylim((0, 20))
plt.tight_layout()
plt.savefig("./img/strong_scaling_v256.png")

# # loading data
# data = np.loadtxt("strong_v1024.dat", delimiter=' ', skiprows=1)
# n2 = data[:,0]
# par2 = data[:,1]
# apar2 = data[:,2]
# # calculating speedup
# par_su2 = par2[0]/par2
# apar_su2 = apar2[0]/apar2
# # plotting weak scaling
# fig,ax = plt.subplots()
# ax.plot(n2, par_su2, label = "Parallelized")
# ax.plot(n2, apar_su2, label = "Par + Async")
# ax.set_title("Speedup, vector size 1024")
# ax.set_xlabel("Number of gangs")
# ax.set_xticks(n2[1:])
# ax.set_ylabel("Speedup")
# ax.legend()
# plt.tight_layout()
# plt.savefig("./img/strong_scaling_v1024.png")

# loading data
data = np.loadtxt("weak_v256.dat", delimiter=' ', skiprows=1)
n_w = data[:,0]
par_w = data[:,1]
apar_w = data[:,2]

# calculating speedup
par_w_su = par_w[0]/par_w
apar_w_su = apar_w[0]/apar_w


# plotting weak scaling
fig,ax = plt.subplots()
ax.plot(n_w, par_w_su, label = "Parallelized")
ax.plot(n_w, apar_w_su, label = "Par + Async")
ax.set_title("Weak scaling")
ax.set_xlabel("Number of gangs")
ax.set_xticks(n_w)
ax.set_ylabel("Efficiency")
ax.legend()
ax.set_ylim((0, 1.1))
plt.tight_layout()
plt.savefig("./img/weak_scaling_v256.png")