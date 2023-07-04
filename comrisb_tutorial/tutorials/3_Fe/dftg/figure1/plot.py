import h5py
import matplotlib.pyplot as plt


e_skn = []
k_dist = []
k_tick_label = []
k_tick_pos = []
file_list = ["../u0j0/lowh/bands_plot.h5",
        "../u5j0.8/lowh/bands_plot.h5"]

for fname in file_list:
    with h5py.File(fname, "r") as f:
        e_skn.append(f["/e_skn"][()])
        k_dist.append(f["/k_dist"][()])
        k_tick_label.append([s.decode("utf-8")  \
                for s in f["/k_tick_label"][()]])
        k_tick_pos.append(f["/k_tick_pos"][()])


fig, ax = plt.subplots(figsize=(4, 4))
for n in range(e_skn[0].shape[2]):
    if n == 0:
        ax.plot(k_dist[0], e_skn[0][0, :, n], 'k--', label="LDA")
        ax.plot(k_dist[1], e_skn[1][0, :, n], 'r-', label="LDA+GRISB")
    else:
        ax.plot(k_dist[0], e_skn[0][0, :, n], 'k--')
        ax.plot(k_dist[1], e_skn[1][0, :, n], 'r-')

ax.axhline(y = 0, ls = ':', lw = 2)
# High-symmetry lines and labels
for x1 in k_tick_pos[0][1:-1]:
    ax.axvline(x = x1, ls = '--')
ax.set_xticks(k_tick_pos[0])
ktick_label = [r"${}$".format(s) for s in k_tick_label[0]]
ax.set_xticklabels(ktick_label)
ax.set_ylabel("E (eV)")
ax.set_xlim(k_dist[0][0], k_dist[0][-1])
ax.set_ylim(-10, 8)
ax.legend()

fig.tight_layout()
plt.show()
fig.savefig("FePMBands.pdf")

