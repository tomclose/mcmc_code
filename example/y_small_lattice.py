import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import matplotlib.pyplot as plt
import y_conjugacy_classes as ycc
import state as st

size = 4

s = st.ToricLattice(size)

@ycc.memoize
def n_to_class(n):
    return ycc.n_to_class(s, n)

@ycc.memoize
def n_to_syndrome(n):
    return ycc.n_to_syndrome(s, n)

def synd_classes():
    synd_classes = {}
    for i in range(2**(size**2)):
        synd = n_to_syndrome(i)
        if (synd,0) in synd_classes:
            # find new class
            c = n_to_class(i^synd_classes[(synd, 0)][0])
            if (synd, c) in synd_classes:
                synd_classes[(synd, c)].append(i)
            else: #new class for that sydrome
                synd_classes[(synd, c)] = [i]
        else: #new syndrome
            synd_classes[(synd, 0)] = [i]
    return synd_classes




