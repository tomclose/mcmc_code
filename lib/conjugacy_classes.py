import numpy as np
import matplotlib.pyplot as plt
import matplotlib as matplotlib # wtf? it should work without this

def bitsum(x):
    return sum(int(i) for i in bin(x)[2:])

def to_n(state):
    # WARNING - only represents z errors
    ans = 0
    for i,j in state.qubit_indices:
        ans = ans*2 + state.array[i,j]
    return ans

class ConjugacyClasses:
    """ Generates a full set of states and categorises
    them into conjugacy classes.

    """
    def __init__(slf, lattice):
        n = slf.n = lattice.L/2

        s = lattice
        print("n=", n, "2**(n**2-1)=", 2**(n**2-1))


        s2 = s.copy()

        stab_gens = slf.stab_gens = [to_n(s.copy_onto(s2).apply_stabiliser(i, j)) for i,j in s.z_stabiliser_indices[:-1]]

        z_hor = slf.z_hor = to_n(s.copy_onto(s2).apply_hor_z())
        z_vert = slf.z_vert = to_n(s.copy_onto(s2).apply_vert_z())
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




