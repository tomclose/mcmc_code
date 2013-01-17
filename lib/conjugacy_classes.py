import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib as matplotlib # wtf? it should work without this

def bitsum(x):
    return sum(int(i) for i in bin(x)[2:])

def to_n(state):
    # WARNING - only represents z errors
    ans = 0
    l = 0
    for i,j in state.qubit_indices:
        if state._array[i,j] > 1:
            raise RuntimeError('to_n passed a[{0}][{1}] = {2}: should only be used for z errors (0 or 1)'.format(i,j,state._array[i,j]))
        ans += (state._array[i,j] << l)
        l += 1
    return ans

def set_state(state, n):
    l = 0
    for i, j in state.qubit_indices:
        if n & (1 << l):
            state.set_qubit(i, j, state.Z)
        else:
            state.set_qubit(i, j, 0)
        l += 1
    return state

class ConjugacyClasses:
    """ Generates a full set of states and categorises
    them into conjugacy classes.

    """
    def __init__(slf, lattice):
        n = slf.n = lattice.L/2

        s = lattice
        print("n=", n, "2**(n**2-1)=", 2**(n**2-1))


        s2 = s.copy()

        stab_gens = slf.stab_gens = [to_n(s.copy_onto(s2).apply_stabiliser(i, j)) for i,j in s.x_stabiliser_indices[:-1]]

        z_hor = slf.z_hor = to_n(s.copy_onto(s2).change_class(lattice.HOR_Z))
        z_vert = slf.z_vert = to_n(s.copy_onto(s2).change_class(lattice.VERT_Z))
        z_hor_vert = slf.z_hor_vert = z_hor ^ z_vert

        print("stab_gens found")

        # builds up stabiliser labelled by n by combining
        # single stabilisers
        def stabiliser(stab_gens, n):
            stab = 0
            bitstr = bin(n)[2:] # change n to binary string
            m = len(bitstr)
            for i in range(m):
                stab ^= int(bitstr[-i-1]) * stab_gens[i]
            return stab

        # complete set of stabilisers
        slf.stabs = stabs = [stabiliser(stab_gens, m) for m in range(2**(n**2-1))]

        print("stabs found")

        # allocate hash to store which orbit each config
        # belongs to 
        slf.record = record = dict([m, None] for m in range(2**(2*n**2)))

        print("dict allocated")

        # table to record the orbtis
        slf.results = results = np.zeros((2**(n**2+1), 2**(n**2-1)), dtype = 'int')

        print("result array allocated")

        m = 0

        for i in range(0, 2**(n**2-1)):
            # find smallest remaining
            while record[m] is not None:
                m += 1
            new_orbit = [x^m for x in stabs]
            new_z_hor = [x^z_hor for x in new_orbit]
            new_z_vert =[x^z_vert for x in new_orbit] 
            new_z_hor_vert =[x^z_hor_vert for x in new_orbit]
            results[4*i, :] = new_orbit
            results[4*i+1, :] = new_z_hor
            results[4*i+2, :] = new_z_vert
            results[4*i+3, :] = new_z_hor_vert
            for j in new_orbit:
                record[j] = 4*i
                record[j^z_hor] = 4*i+1
                record[j^z_vert] = 4*i+2
                record[j^z_hor_vert] = 4*i+3

            if i%2**10==0:
                print(i, 2**(n**2-1))

        b = np.vectorize(bitsum)

        x = b(results)

        for r in x:
            r[:] = np.sort(r)[:]

        slf.hist = np.array([np.histogram(a, bins=np.arange(-0.5, 2*n**2+1.5, 1))[0] for a in x])

    @staticmethod
    def imshowx(x):
        y_max, x_max = np.shape(x)
        colors = [('white')] + [(plt.cm.jet(i)) for i in xrange(1,256)]
        new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)
        im = plt.imshow(x, cmap=new_map, interpolation='nearest', aspect='auto', extent=(0,x_max, 0, y_max) )
        ax = plt.gca()
        ax.yaxis.set_minor_locator(plt.MultipleLocator(4))
        plt.grid(which='minor')
        return im

    def show(s):
        b = np.vectorize(bitsum)
        s.imshowx(b(s.results))

    def show_hist(s):
        s.imshowx(s.hist)


    def lookup(s, n):
        # returns the conjugacy classes for s, X_vS, X_hS, X_vX_hS
        i = s.record[n]
        v = s.record[s.z_vert^n]
        h = s.record[s.z_hor^n]
        vh = s.record[s.z_hor_vert^n]
        return s.results[[i, v, h, vh]]


    def lookup_hist(s, n):
        # returns the conjugacy classes for s, X_vS, X_hS, X_vX_hS
        i = s.record[n]
        v = s.record[s.z_vert^n]
        h = s.record[s.z_hor^n]
        vh = s.record[s.z_hor_vert^n]
        return s.hist[[i, v, h, vh]]


    @staticmethod
    def hist_row_prob(hist_row, p):
        q = p/(1-p)
        return sum(n*q**i for (i, n) in enumerate(hist_row))

    #def sort_errors(s, p):
        #for i in range(len(s.hist)/4):
            #r = s.hist[4*i:4*(i+1)]

            #r[:] = r[np.argsort([-1*s.hist_row_prob(row, p) for row in r])]
        #return s.hist

    def total_prob(s, p):
        n_qubits = (s.n/2)**2
        h = s.sort_errors(s.hist, p)
        return (1-p)**n_qubits * sum(s.hist_row_prob(y, p) for y in h[::4]) 

    @staticmethod
    def hist_row_av_n(hist_row, p):
        q = p/(1-p)
        e = sum(i*n*q**i for (i, n) in enumerate(hist_row))
        totp = sum(n*q**i for (i, n) in enumerate(hist_row))
        return e*1.0/totp

        

# the probability of the correct identification is the probability it is the top group x the probability you 
# identify it as the top group. The probability you identify it as the top group is roughly proportional 

# in fact assume we can identify the best group arbitrarily (we could decide the result for each syndrome beforehand)
# then the probability is just the sum of the probabilities of the best group out of each

def stab_gen_effect(state, i, j):
    state.apply_stabiliser(i, j)
    n = to_n(state)
    state.apply_stabiliser(i, j) # undo
    return n

def stab_gens(s):
    if s.n_errors() > 0:
        raise RuntimeError("Must be passed a clean state")
    return [stab_gen_effect(s, i, j) for i, j in s.z_stabiliser_indices]

def stab_effect(stab_gens, n):
    """ Returns the stabiliser effect for stabiliser
    described by n
    NB: 0 <= n < 2**m, but effect will only be unique
        for 0 <= n < 2**(m-1)
    """
    stab = 0
    bitstr = bin(n)[2:] # change n to binary string
    m = len(bitstr)
    for i in range(m):
        stab ^= int(bitstr[m-1-i]) * stab_gens[i]
    return stab

def stabilisers(s):
    gens = stab_gens(s)
    n = len(gens)
    return [stab_effect(gens, x) for x in range(2**(n - 1))]

def orbit(stabs, n):
    return [n^s for s in stabs]

def conj_class_gens(s):
    vert_z = to_n(s.change_class(s.VERT_Z))
    s.change_class(s.VERT_Z)
    hor_z = to_n(s.change_class(s.HOR_Z))
    hor_vert_z = to_n(s.change_class(s.VERT_Z))
    s.change_class(s.VERT_Z + s.HOR_Z) # reset
    return [0, vert_z, hor_z, hor_vert_z]

def matching(syndrome_n, state):
    stab_ind = state.x_stabiliser_indices #=> (1,1), (1,3), (3,1), (3,3)
    n_stabs = len(stab_ind) #=> 4
    synd = []
    # map the syndrome_n onto the first n_stab -1 indices
    for i in range(n_stabs - 1): # eg 0, 1, 2
        if syndrome_n & ( 1 << (n_stabs - 1 - i)): # eg 
            synd.append(stab_ind[i])
    if bitsum(syndrome_n) & 1: # odd number of stabilisers
        synd.append(stab_ind[-1])
    match = state.__class__.from_syndrome(state.L, synd)
    return to_n(match)

def error_dist(orbit):
    dist = {}
    for x in orbit:
        b = bitsum(x)
        try:
            dist[b] += 1
        except KeyError: # first time we've seen it
            dist[b] = 1
    return dist

def hist_row(error_dist, max_poss):
    return [error_dist.get(i, 0) for i in range(max_poss + 1)]


def state_list(s):
    stabs = stabilisers(s)
    conj_classes = conj_class_gens(s)
    for i, stab in enumerate(stabs):
        match = matching(i, s)
        for conj_class in conj_classes:
            yield orbit(stabs, match^conj_class)

def hist_file_rows(s):
    stabs = stabilisers(s)
    conj_classes = conj_class_gens(s)
    for i, stab in enumerate(stabs):
        match = matching(i, s)
        for conj_class in conj_classes:
            orb = orbit(stabs, match^conj_class)
            yield [ stab, conj_class] + hist_row(error_dist(orb), s.L**2/2)

def write_hist_file(s, filename = None):
    if filename is None:
        filename = './data/hist_{}.csv'.format(s.L)
    with open(filename, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for histrow in hist_file_rows(s):
            writer.writerow(histrow)

def read_hist_file_rows(filename):
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            yield row

def array_from_file_rows(rows):
    return np.array([r for r in rows], dtype='int')

def class_probabilities(row_array, p):
    # row_array has format [syndrome, logical_error, count(n_errors = 0), count(n_errors = 1), .... ]
    num_classes, total_errors = np.shape(row_array[:, 2:])
    if p == 0:
        q = 0
    else:
        q = p/(1-p)
    prob_row = q ** np.arange(total_errors)
    class_probs = (1-p)**(total_errors-1) * np.dot(row_array[:, 2:], prob_row)
    results = np.zeros((num_classes, 2))
    results[:, 0] = row_array[:, 0]
    results[:, 1] = class_probs
    return results

