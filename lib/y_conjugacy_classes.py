import numpy as np
import state as st

#http://wiki.python.org/moin/PythonDecoratorLibrary#Memoize
class memoize(dict):
    def __init__(self, func):
        self.func = func

    def __call__(self, *args):
        return self[args]

    def __missing__(self, key):
        result = self[key] = self.func(*key)
        return result

@memoize
def bitsum(x):
    """Return the number of 1s in the binary representation of x"""
    return bin(x).count('1')

def set_state_errors(state, n):
    """Returns state with error configuration represented by n"""
    l = 0
    for i, j in state.qubit_indices:
        if n & (1 << l):
            state.set_qubit(i, j, state.Y)
        else:
            state.set_qubit(i, j, 0)
        l += 1
    return state

def errors_to_n(state):
    """
      Maps the state's qubit error configuration to an integer, n.

      The qubit in the top left corresponds to the least significant bit in n.

    """
    ans = 0
    l = 0
    for i,j in state.qubit_indices:
        error = state._array[i, j]
        if error == state.X:
            raise RuntimeError('to_n passed a[{0}][{1}] = {2}: should only be used for z errors (0 or 1)'.format(i,j,state._array[i,j]))
        elif error == state.Z:
            raise RuntimeError('to_n passed a[{0}][{1}] = {2}: should only be used for z errors (0 or 1)'.format(i,j,state._array[i,j]))
        elif error == state.Y:
            ans += (1 << l)
        else:
            ans += (0 << l)
        l += 1
    return ans

def syndrome_to_n(state):
    """
      Maps the state's syndrome an integer, n.

      The stabiliser in the top left corresponds to the least significant bit in n.

    """
    ans = 0
    l = 0
    for i,j in state.stabiliser_indices:
        stab_val = state._array[i, j]
        if stab_val == 0:
            ans += (0 << l)
        elif stab_val == 1:
            ans += (1 << l)
        else:
            raise RuntimeError('syndrome_to_n passed a[{0}][{1}] = {2}: should only be 0 or 1'.format(i,j,state._array[i,j]))
        l += 1
    return ans

def n_to_class(state, n):
    """Returns class for n (assuming it is in the code space).

       State is a dummy state of the correct size."""
    set_state_errors(state, n)
    return st.ToricLattice.compare(state)

def n_to_syndrome(state, n):
    set_state_errors(state, n)
    state.syndrome(refresh=True)
    return syndrome_to_n(state)

def error_dist(orbit):
    dist = {}
    for x in orbit:
        b = bitsum(x)
        try:
            dist[b] += 1
        except KeyError: # first time we've seen it
            dist[b] = 1
    return dist

def hist_row(orbit, max_poss):
    dist = error_dist(orbit)
    return [dist.get(i, 0) for i in range(max_poss +1)]


def hist_array(synd_clss, size):
    max_poss = size**2/2
    array = [[k[0], k[1]] + hist_row(v, max_poss) for k, v in synd_clss.iteritems()]
    array = np.array(array)
    # collect different classes together
    array = array[np.argsort(array[:, 0])]
    return array

def class_probabilities(row_array, p):
    # row_array has format [syndrome, logical_error, count(n_errors = 0), count(n_errors = 1), .... ]
    num_classes, total_errors = np.shape(row_array[:, 2:])
    if p < 0.5: # rules out p = 1
        q = p/(1-p)
        prob_row = (1-p)**(total_errors-1) * q ** np.arange(total_errors)
    else: # p != 0
        q = (1-p)/p
        prob_row = p **(total_errors-1) * q ** np.arange(total_errors-1, -1, -1)
    class_probs = np.dot(row_array[:, 2:], prob_row)
    results = np.zeros((num_classes, 2))
    results[:, 0] = row_array[:, 0]
    results[:, 1] = class_probs
    return results

# memoize
def max_class_probs(row_array, p):
    probs = class_probabilities(row_array, p)
    x = probs[:, 1].reshape((-1, 4))
    return np.max(x, axis = 1)

def success_probability(row_array, p):
    return np.sum(max_class_probs(row_array, p))

