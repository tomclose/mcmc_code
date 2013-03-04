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


