import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import conjugacy_classes as cc
import pickle

#http://wiki.python.org/moin/PythonDecoratorLibrary#Memoize
class memoize(dict):
    def __init__(self, func):
        self.func = func

    def __call__(self, *args):
        return self[args]

    def __missing__(self, key):
        result = self[key] = self.func(*key)
        return result



rows2 = cc.read_hist_file_rows('./data/hist_2.csv')
rows4 = cc.read_hist_file_rows('./data/hist_4.csv')
rows6 = cc.read_hist_file_rows('./data/hist_6.csv')
rows8 = cc.read_hist_file_rows('./data/hist_8.csv')

a2 = cc.array_from_file_rows(rows2)
a4 = cc.array_from_file_rows(rows4)
a6 = cc.array_from_file_rows(rows6)
a8 = cc.array_from_file_rows(rows8)

def p2n(qubit_error_p, stab_error_p):
    return cc.noisy_prob(a2, 2, qubit_error_p, stab_error_p)

def p4n(qubit_error_p, stab_error_p):
    return cc.noisy_prob(a4, 4, qubit_error_p, stab_error_p)

@memoize
def p6n(qubit_error_p, stab_error_p):
    return cc.noisy_prob(a6, 6, qubit_error_p, stab_error_p)

# load in memoized values
with open('./data/p6n.pickle', 'rb') as f:
    vals = pickle.load(f)
    p6n.update(vals)


#def p8n(qubit_error_p, stab_error_p):
    #return cc.noisy_prob(a8, 8, qubit_error_p, stab_error_p)

qq = pp = np.linspace(0, 0.1, 51)

PP, QQ = np.meshgrid(pp, qq)

Z2 = np.array([[p2n(p, q) for p in pp] for q in qq])
Z4 = np.array([[p4n(p, q) for p in pp] for q in qq])
Z6 = np.array([[p6n(p, q) for p in pp] for q in qq])
#Z8 = np.array([[p8(p, q) for p in pp] for q in qq])

bare_qubit = PP

#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)

fontsize=16

ax.contour(QQ, PP, Z6 - (1 - PP), 100)
ax.contour(QQ, PP, Z6 - (1 - PP), [0], colors = ('black'), linewidths=(3))
plt.xlabel('Probability that stabiliser lies', fontsize=fontsize)
plt.ylabel('Qubit error probability', fontsize=fontsize)
plt.show()
plt.savefig('x_lying.pdf', format='pdf')
plt.savefig('../writeup/assets/x_lying.pdf', format='pdf')
