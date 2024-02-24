import numpy as np
import matplotlib.pyplot as plt

# loading data
data = np.loadtxt("perf.dat", delimiter=' ', skiprows=1, max_rows=6)
n = data[:,0]
total_time = data[:,1]
task_time = data[:,2]

# calculating speedup
total_su = total_time[0]/total_time
task_su = task_time[0]/task_time

# plotting absolute performance
fig,ax = plt.subplots()
ax.plot(n, task_time)
ax.set_title("Time per task")
ax.set_xlabel("Number of processes")
ax.set_xticks(n)
ax.set_ylabel("Microseconds")
plt.tight_layout()
plt.savefig("./img/task_perf.png")

# plotting relative performance
fig,ax = plt.subplots()
ax.plot(n, task_su)
ax.set_title("Speedup")
ax.set_xlabel("Number of processes")
ax.set_xticks(n)
ax.set_ylabel("Speedup")
plt.tight_layout()
plt.savefig("./img/task_speedup.png")


### loading data from multiple compute nodes
# data_all = np.loadtxt("perf.dat", delimiter=' ', skiprows=1, max_rows=8)
# n_all = data_all[:,0]
# total_time_all = data_all[:,1]
# task_time_all = data_all[:,2]
# task_su_all = task_time_all[0]/task_time_all

# # plotting relative performance
# fig,ax = plt.subplots()
# ax.plot(n_all, task_su_all)
# ax.set_title("Speedup")
# ax.set_xlabel("Number of processes")
# ax.set_xticks(n_all)
# ax.set_ylabel("Speedup")
# plt.tight_layout()
# plt.savefig("./img/task_speedup_all.png")


### calculating parallel fraction
workers = n[1:] - 1
p_frac = (-1 + 1/total_su[1:])/(-1 + 1/workers)
print("Workers:", workers)
print("Speedup:", np.around(total_su[1:], decimals=2))
print("Parallel Frac:", np.around(p_frac[1:], decimals=3))

### theoretical scaling curve
def amdahl(n, p):
    return 1 / ((1-p) + p/n)

workers_theo = np.full(12, 2) ** np.arange(0, 12)
speedup_theo = amdahl(workers_theo, 0.993)
#speedup_theo2 = amdahl(workers_theo, 0.988)

fig,ax = plt.subplots()
ax.plot(workers_theo, speedup_theo, label = "0.7% serial")
#ax.plot(workers_theo, speedup_theo2, label = "1.2% serial")
ax.set_title("Theoretical speedup")
ax.set_xlabel("Number of processes")
ax.set_xticks([16, 128, 256, 512, 1028, 2048])
ax.legend()
ax.set_ylabel("Speedup")

plt.tight_layout()
plt.savefig("./img/theo_speedup.png")