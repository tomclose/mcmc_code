import numpy as np
import matplotlib.pyplot as plt

class State:
    def __init__(self, L, p):
        self.L = L
        self.array = np.zeros((L,L), dtype=('bool', 1))
        self.p_val = p

    def show(self):
        fig = plt.gcf()# or plt.figure()
        fig.clf()
        ax = fig.add_subplot(111)
        # Decide how to show array:
        #  - 0 is a non-firing stabiliser
        #  - 1 is a non-error qubit
        #  - 2 is a firing stabiliser
        #  - 3 is an error qubit
        show_array = np.zeros((self.L, self.L), dtype='int')
        show_array += 2* self.array
        show_array += self.qubit_mask()
        cax = ax.imshow(show_array, interpolation='nearest', vmax=3)
        cbar = plt.colorbar(cax, ticks=[0, 1, 2, 3])
        cbar.ax.set_yticklabels(['stab', 'qubit', 'firing stab', 'error qubit'])
        plt.show() # in case not already drawn
        plt.draw() # refresh if drawn
        return show_array

    # Array:
    #   - If sum is odd then it represents a qubit (X_error?, Y_error?)
    #   - If sum is even then
    #       - If even + even it represents an X stabiliser
    #       - If odd + odd it represents a Z stabiliser
    # Actually - do the simple case first if sum is odd just Z_error?

    def qubit_mask(self):
        a = np.zeros((self.L, self.L), dtype='bool')
        for i in range(self.L):
            for j in range(self.L):
                if (j+i)%2 == 1:
                    a[i, j] = 1
        return a

    def x_stabiliser_mask(self):
        """Returns a boolean array of x stabiliser positions.

        The x stabilisers occur on (i,j) where both i and j are even
        """
        a = np.zeros((self.L, self.L), dtype='bool')
        for i in range(0,self.L, 2):
            for j in range(0, self.L, 2):
                    a[i, j] = 1
        return a

    def z_stabiliser_mask(self):
        """Returns a boolean array of z stabiliser positions.

        The z stabilisers occur on (i,j) where both i and j are odd
        """
        a = np.zeros((self.L, self.L), dtype='bool')
        for i in range(1,self.L, 2):
            for j in range(1, self.L, 2):
                    a[i, j] = 1
        return a
    
    def generate_x_syndrome(self):
        for i in range(0,self.L, 2):
            for j in range(0, self.L, 2):
                a = self.array
                nl, nr, nt, nb = a[i-1, j], a[(i+1)%self.L, j], a[i, j-1], a[i, (j+1)%self.L]
                n = 0 + nl + nr + nt + nb # the zero enforces upcasting to int
                ans = n % 2
                if ans == 0:
                    #print("True", (i,j), [nl, nr, nt, nb], n, ans)
                    a[i, j] = 0
                else:
                    #print("False", (i,j), [nl, nr, nt, nb], n,ans)
                    a[i, j] = 1
        return a

    def likelihood(self):
        n = np.sum(self.array[self.qubit_mask()])
        N = np.sum(self.qubit_mask())
        p = self.p_val
        return p**n * (1-p)**(N-n)

    def generate_errors(self):
        n_qubits = np.sum(self.qubit_mask())
        errors = np.random.rand(n_qubits) < self.p_val
        self.array[self.qubit_mask()] = errors

    def generate_next(self):
        # pick a random z site
        xs, ys = np.where(self.z_stabiliser_mask()==1)
        n = len(xs)
        r = np.random.random_integers(0, n-1)
        x, y = xs[r], ys[r]
        # apply the stabilizer at that point
        a = self.array
        L = np.shape(a)[0]
        a[x-1, y] = not a[x-1, y] # logical flip
        a[x, y-1] = not a[x, y-1] # logical flip
        a[(x+1)%L, y] = not a[(x+1)%L, y] # logical flip
        a[x, (y+1)%L] = not a[x, (y+1)%L] # logical flip


