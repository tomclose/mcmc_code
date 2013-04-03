import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import matplotlib.pyplot as plt
import y_conjugacy_classes as ycc
import state as st

s2 = st.ToricLattice(2)
sc2 = ycc.synd_classes(s2)
h2 = ycc.hist_array(sc2, 2)

s4 = st.ToricLattice(4)
sc4 = ycc.synd_classes(s4)
h4 = ycc.hist_array(sc4, 4)

s6 = st.ToricLattice(6)
sc6 = ycc.synd_classes(s6)
h6 = ycc.hist_array(sc6, 6)

def p2(p):
    return ycc.success_probability(h2, p)
def p4(p):
    return ycc.success_probability(h4, p)
def p6(p):
    return ycc.success_probability(h6, p)

pp = np.linspace(0, 1, 101)

fontsize=16

plt.clf()
#plt.plot(pp, 1 - pp, label='bare qubit')
#plt.plot(pp, (1 - pp)**2, label='2 bare qubits')
plt.plot(pp, [p2(p) for p in pp], label='2-code')
plt.plot(pp, [p4(p) for p in pp], label='4-code')
plt.plot(pp, [p6(p) for p in pp], label='6-code')
plt.legend()
plt.xlabel('$p$', fontsize=fontsize)
plt.ylabel('$P_d$', fontsize=fontsize)
plt.show()
plt.savefig('y_truthful.pdf', format='pdf')
plt.savefig('../writeup/assets/y_truthful.pdf', format='pdf')
