import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import matplotlib.pyplot as plt
import conjugacy_classes as cc

rows2 = cc.read_hist_file_rows('./data/hist_2.csv')
rows4 = cc.read_hist_file_rows('./data/hist_4.csv')
rows6 = cc.read_hist_file_rows('./data/hist_6.csv')
rows8 = cc.read_hist_file_rows('./data/hist_8.csv')

a2 = cc.array_from_file_rows(rows2)
a4 = cc.array_from_file_rows(rows4)
a6 = cc.array_from_file_rows(rows6)
a8 = cc.array_from_file_rows(rows8)

def p2(p):
    return cc.success_probability(a2, p)

def p4(p):
    return cc.success_probability(a4, p)

def p6(p):
    return cc.success_probability(a6, p)

def p8(p):
    return cc.success_probability(a8, p)

pp = np.linspace(0, 1, 101)

fontsize=16

plt.clf()
plt.plot(pp, 1 - pp, label='bare qubit')
plt.plot(pp, (1 - pp)**2, label='2 bare qubits')
plt.plot(pp, [p2(p) for p in pp], label='2-code')
plt.plot(pp, [p4(p) for p in pp], label='4-code')
plt.plot(pp, [p6(p) for p in pp], label='6-code')
plt.plot(pp, [p8(p) for p in pp], label='8-code')
plt.legend()
plt.xlabel('Qubit error probability', fontsize=fontsize)
plt.ylabel('Decoding success probability', fontsize=fontsize)
plt.show()
plt.savefig('x_truthful.pdf', format='pdf')
plt.savefig('../writeup/assets/x_truthful.pdf', format='pdf')
