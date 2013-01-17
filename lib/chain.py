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
        #p = current_state.likelihood()
        #q = new.likelihood()
        rel_prob = current_state.relative_prob(new) # q/p
        #print(rel_prob)
        #print(p, q)
        # r = q/p not robust, instead of X < r do p*X < q
        if np.random.rand() < rel_prob: # includes the r>1 case
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

def in_jumps_of(n, c):
    return it.islice(c, n-1, None, n)

def first_n(n, c):
    return it.islice(c, 0, n)

def path_set(*paths):
    return it.izip(*paths)

def totals(p):
    n_changes = 0
    count = 0
    total_n = 0
    for current, next_s, changed in p:
        n_changes += changed
        count+=1
        total_n += current.n_errors()
        yield (total_n, n_changes, count)


def average_err(p):
    n_changes = 0
    count = 0
    total_n = 0
    for current, next_s, changed in p:
        n_changes += changed
        count+=1
        total_n += current.n_errors()
        yield (total_n, n_changes, count)
