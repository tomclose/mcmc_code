import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import y_conjugacy_classes as ycc
import state as st

s4 = st.ToricLattice(4)
sc4 = ycc.synd_classes(s4)
h4 = ycc.hist_array(sc4, 4)

def p4n(p, q):
    return ycc.small_noisy_prob(h4, 4, p, q)

pp = np.linspace(0, 0.2, 81)
qq = np.linspace(0, 0.12, 49)

PP, QQ = np.meshgrid(pp, qq)

Z4 = np.array([[p4n(p, q) for p in pp] for q in qq])

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)

fontsize=16

Z = Z4 - (1 - PP)**2
cont =  ax.contour(QQ, PP, Z, [0.06, 0.05, 0.04, 0.03, 0.02, 0.01], colors=('black'))
plt.clabel(cont, inline=1, fontsize=fontsize)
cont2 = ax.contour(QQ, PP, Z, [0], colors = ('black'), linewidths=(3))
plt.clabel(cont2, inline=1, fontsize=fontsize)
plt.xlabel('$q$', fontsize=fontsize)
plt.ylabel('$p$', fontsize=fontsize)
plt.show()
plt.savefig('y_lying.pdf', format='pdf')
plt.savefig('../writeup/assets/y_lying.pdf', format='pdf')
