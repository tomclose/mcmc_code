import sys
import os
# this path stuff is so retarded
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import state
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

n = 3
print("n=", n, "2**(n**2-1)=", 2**(n**2-1))


s = state.State(2*n, 0.3)
s2 = s.copy()

stab_gens = [s.copy_onto(s2).apply_stabiliser(i, j).to_n() for i,j in s.x_stabiliser_indices()[:-1]]

z_hor = s.copy_onto(s2).generate_horizontal_z_error().to_n()
z_vert = s.copy_onto(s2).generate_vertical_z_error().to_n()
z_hor_vert = z_hor ^ z_vert

print("stab_gens found")

def stabiliser(stab_gens, n):
    stab = 0
    bitstr = bin(n)[2:]
    m = len(bitstr)
    for i in range(m):
        stab ^= int(bitstr[-i-1]) * stab_gens[i]
    return stab

stabs = [stabiliser(stab_gens, m) for m in range(2**(n**2-1))]

print("stabs found")

record = dict([m, None] for m in range(2**(2*n**2)))

print("dict allocated")

results = np.zeros((2**(n**2+1), 2**(n**2-1)), dtype = 'int')

print("result array allocated")

#for x in stabs:
    #record[x] = 0

#results[0, :] = stabs[:]
#results[1, :] = stabs[:]^z_hor
#results[1, :] = stabs[:]^z_hor
#results[1, :] = stabs[:]^z_hor

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

b = np.vectorize(state.bitsum)

x = b(results)

for r in x:
    r[:] = np.sort(r)[:]

def imshowx(x):
    y_max, x_max = np.shape(x)
    colors = [('white')] + [(plt.cm.jet(i)) for i in xrange(1,256)]
    new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)
    im = plt.imshow(x, cmap=new_map, interpolation='nearest', aspect='auto', extent=(0,x_max, 0, y_max) )
    ax = plt.gca()
    ax.yaxis.set_minor_locator(plt.MultipleLocator(4))
    plt.grid(which='minor')
    return im


h = np.array([np.histogram(a, bins=np.arange(-0.5, 2*n**2+1.5, 1))[0] for a in x])

def hist_row_prob(hist_row, p):
    q = p/(1-p)
    return sum(n*q**i for (i, n) in enumerate(hist_row))

def sort_errors(hist, p):
    for i in range(len(hist)/4):
        r = hist[4*i:4*(i+1)]

        r[:] = r[np.argsort([-1*hist_row_prob(row, p) for row in r])]
    return hist

def total_prob(hist, p, n_qubits):
    h = sort_errors(hist, p)
    return (1-p)**n_qubits * sum(hist_row_prob(y, p) for y in h[::4]) 


# the probability of the correct identification is the probability it is the top group x the probability you 
# identify it as the top group. The probability you identify it as the top group is roughly proportional 

# in fact assume we can identify the best group arbitrarily (we could decide the result for each syndrome beforehand)
# then the probability is just the sum of the probabilities of the best group out of each




