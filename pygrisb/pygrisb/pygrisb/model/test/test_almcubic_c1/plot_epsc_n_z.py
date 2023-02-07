import numpy
import matplotlib.pyplot as plt


data = numpy.loadtxt("ref_epsc_n_z.dat").T

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('$\epsilon_{c}$')
ax1.set_ylabel('$n_{e}$', color=color)
ax1.plot(data[0], data[1], color=color)
ax1.tick_params(axis='y', labelcolor=color)
# instantiate a second axes that shares the same x-axis
ax2 = ax1.twinx()

color = 'tab:blue'
ax2.set_ylabel('$Z_{small}$', color=color)
ax2.plot(data[0], data[2], color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.tight_layout()
plt.show()
