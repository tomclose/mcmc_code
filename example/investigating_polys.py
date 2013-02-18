import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import conjugacy_classes as cc
import matplotlib.pyplot as plt
from collections import defaultdict

rows2 = cc.read_hist_file_rows('./data/hist_2.csv')
rows4 = cc.read_hist_file_rows('./data/hist_4.csv')
rows6 = cc.read_hist_file_rows('./data/hist_6.csv')
rows8 = cc.read_hist_file_rows('./data/hist_8.csv')

a2 = cc.array_from_file_rows(rows2)
a4 = cc.array_from_file_rows(rows4)
a6 = cc.array_from_file_rows(rows6)
a8 = cc.array_from_file_rows(rows8)

d2 = defaultdict(int)
for e in a2:
    d2[tuple(e[2:])]+=1
d4 = defaultdict(int)
for e in a4:
    d4[tuple(e[2:])]+=1
d6 = defaultdict(int)
for e in a6:
    d6[tuple(e[2:])]+=1
d8 = defaultdict(int)
for e in a8:
    d8[tuple(e[2:])]+=1



b = a4[:,2:]

def plot(n):
    filename = './data/hist_{}.csv'.format(n)
    rows = cc.read_hist_file_rows(filename)
    a = cc.array_from_file_rows(rows)
    b = a[:, 2:] # get rid of first (identity) row

    def chi(x):
        xx = [x**i for i in range(2*(n/2)**2+1)] 
        return np.dot(b, xx)

    xx = np.linspace(0, 0.4, 51)

    def chi_vals(xx):
        return np.array([chi(x) for x in xx])

    chivals = chi_vals(xx)

    for i in range(1, 2*((n/2)**2+1)):
        plt.plot(xx, chivals[:, i])


