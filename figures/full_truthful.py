import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import matplotlib.pyplot as plt
import full_conj_classes as fcc

a2 = fcc.read_hist_file(2)
a4 = fcc.read_hist_file(4)
a6 = fcc.read_hist_file(6)

def p2(p):
    return fcc.success_probability(a2, p)

def p4(p):
    return fcc.success_probability(a4, p)

def p6(p):
    return fcc.success_probability(a6, p)

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
plt.ylabel('P_d', fontsize=fontsize)
plt.show()
plt.savefig('full_truthful.pdf', format='pdf')
plt.savefig('../writeup/assets/full_truthful.pdf', format='pdf')
