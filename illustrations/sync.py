import plotext as plt
import numpy as np

p = 1.2

gs = np.logspace(1, 5)
dn_dg = (gs / gs[0]) ** (-p) * np.exp(-(gs / 1e3))

plt.theme("pro")
plt.plot_size(80, 40)
plt.plot(np.log10(gs), np.log10(dn_dg))
plt.xlim(0.5, 4)
plt.ylim(-7, 1)
plt.show()
