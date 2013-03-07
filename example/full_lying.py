import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import full_conj_classes as fcc
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

def save_memoized(fcn, filename):
    path = './data/' + filename + '.csv'
    with open(path, 'wb') as f:
        writer = csv.writer(f)
        for (k1, k2), v in fcn.items():
            writer.writerow([k1, k2, v])

def load_memoized(fcn, filename):
    path = './data/' + filename + '.csv'
    with open(path, 'rb') as f:
        reader = csv.reader(f)
        for k1, k2, v in reader:
            fcn[(float(k1), float(k2))] = float(v)


a2 = fcc.read_hist_file(2)
a4 = fcc.read_hist_file(4)
a6 = fcc.read_hist_file(6)


def p2n(qubit_error_p, stab_error_p):
    return fcc.noisy_prob(a2, 2, qubit_error_p, stab_error_p, 2)

def p4n(qubit_error_p, stab_error_p):
    return fcc.noisy_prob(a4, 4, qubit_error_p, stab_error_p, 8)

@memoize
def p6n(qubit_error_p, stab_error_p):
    return fcc.noisy_prob(a6, 6, qubit_error_p, stab_error_p, 4)

@memoize
def p6n2(qubit_error_p, stab_error_p):
    return fcc.noisy_prob(a6, 6, qubit_error_p, stab_error_p, 2)

load_memoized(p6n2, 'p6n2')


#def p8n(qubit_error_p, stab_error_p):
    #return cc.noisy_prob(a8, 8, qubit_error_p, stab_error_p)

qq = pp = np.linspace(0, 0.2, 41)

PP, QQ = np.meshgrid(pp, qq)

Z2 = np.array([[p2n(p, q) for p in pp] for q in qq])
Z4 = np.array([[p4n(p, q) for p in pp] for q in qq])
Z6 = np.array([[p6n2(p, q) for p in pp] for q in qq])

save_memoized(p6n2, 'p6n2')

#Z8 = np.array([[p8(p, q) for p in pp] for q in qq])




bare_qubit = PP

#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)


ax.contour(QQ, PP, Z6 - (1 - PP), 100)
ax.contour(QQ, PP, Z6 - (1 - PP), [0], colors = ('black'), linewidths=(3))
plt.xlabel('Lying stabiliser')
plt.ylabel('Qubit error')
