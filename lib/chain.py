import numpy as np
import itertools as it
#import state

def path(s0):
    current_state = s0
    next_state = s0.copy()
    changed = False # did current_state change last time

    while(1):
        new = current_state.copy_onto(next_state)
        new.generate_next()
        p = current_state.likelihood()
        q = new.likelihood()
        #print(p, q)
        # r = q/p not robust, instead of X < r do p*X < q
        if p * np.random.rand() < q: # includes the r>1 case
            next_state.copy_onto(current_state)
            changed = True
        else:
            changed = False
        yield (current_state,new,  changed)


def aggregate(c):
    p = 0
    i = 0
    for (s, n, changed) in c:
        i+=1
        p += s.likelihood()
        yield (i, p)

def n_jumps(c, n):
    return it.islice(c, n, None, n)

def first_n(c, n):
    return it.islice(c, 0, n)


def path_set(*paths):
    return it.izip(*paths)


