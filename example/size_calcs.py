from scipy.misc import comb
import numpy as np
import matplotlib.pyplot as plt

a = np.array([[comb(2*L**2, r) for r in range(20)] for L in range(0, 20)])
s = np.cumsum(a, axis=1)

cs = plt.contour(s, [10, 100, 1000, 10000, 100000, 1000000, 10000000])
plt.clabel(cs, inline=1, fontsize=10)

plt.show()

