import numpy as np
import itertools as it
#import state

class path:
    def __init__(self, state):
        self.current_state = state
        self.next_state = state.copy()
        self.changed = False # did current_state change last time

    def __iter__(self):
        return self

    def next(self):
        current_state = self.current_state
        new = current_state.copy_onto(self.next_state)
        new.generate_next()
        # = current_state.likelihood()
        #q = new.likelihood()
        rel_prob = current_state.relative_prob(new) # q/p
        #print(rel_prob)
        #print(p, q)
        # r = q/p not robust, instead of X < r do p*X < q
        if np.random.rand() < rel_prob: # includes the r>1 case
            new.copy_onto(current_state)
            self.changed = True
        else:
            self.changed = False
        return current_state

class modified_path:
    def in_jumps_of(self, n):


def in_jumps_of(n, c):
    return it.islice(c, n, None, n)

def first_n(c, n):
    return it.islice(c, 0, n)

def path_set(*paths):
    return it.izip(*paths)


def average_err(p):
    n_changes = 0
    count = 0
    total_n = 0
    for current, next_s, changed in p:
        n_changes += changed
        count+=1
        total_n += current.n_errors()
        yield (total_n, n_changes, count)
