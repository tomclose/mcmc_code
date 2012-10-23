import sys
import os
# this path stuff is so retarded
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import state
import numpy as np
import matplotlib as plt

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

plt.imshow(x, extent = (0, 2**(n**2-1), 0, 2**(n**2+1)), interpolation='nearest')
ax = plt.gca()

