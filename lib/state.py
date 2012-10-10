import numpy as np
import matplotlib.pyplot as plt

class State:
    """ State is reponsible for holding the configuration
    of a toric lattice - both errors, and stabiliser states.
    
    The lattice has side length L, which is taken to be even
    (so the toric edges match).

    The qubits lie in positions (i,j), where i + j is odd.
    Their state is encoded as:
        - 0: no error
        - 1: a Z flip
        - 2: an X flip
        - 3: a Y flip (X and Z)

    The X plaquettes lie at (i, j), where both i and j are even.
    They fire if surrounded by an odd number of {Z,Y} flips

    The Z plaquettes lie at (i, j), where both i and j are odd.
    The fire if surrounded by and odd numper of {X, Y} flips

    Currently we only deal with Z flips, detected by X plaquettes.

    """
    def __init__(self, L, p):
        self.L = L
        self.array = np.zeros((L,L), dtype='uint8')
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

    def x_stabiliser_indices(self):
        return [(i,j) for i in range(0, self.L, 2) for j in range(0, self.L, 2)]
    def z_stabiliser_indices(self):
        return [(i,j) for i in range(1, self.L, 2) for j in range(1, self.L, 2)]
    def qubit_indices(self):
        return [(i,j) for j in range(0, self.L) for i in range((j+1)%2, self.L, 2)]
    def neighbours(self, i, j):
        return [(i-1, j), (i, j-1), ((i+1)%self.L, j), (i, (j+1)%self.L)]

    def qubit_mask(self):
        a = np.zeros((self.L, self.L), dtype='bool')
        for i,j in self.qubit_indices():
            a[i, j] = 1
        return a

    def x_stabiliser_mask(self):
        """Returns a boolean array of x stabiliser positions.

        The x stabilisers occur on (i,j) where both i and j are even
        """
        a = np.zeros((self.L, self.L), dtype='bool')
        for i,j in self.x_stabiliser_indices():
            a[i, j] = 1
        return a

    def z_stabiliser_mask(self):
        """Returns a boolean array of z stabiliser positions.

        The z stabilisers occur on (i,j) where both i and j are odd
        """
        a = np.zeros((self.L, self.L), dtype='bool')
        for i,j in self.z_stabiliser_indices():
            a[i, j] = 1
        return a
    
    def generate_x_syndrome(self):
        coords = []
        for i, j in self.x_stabiliser_indices():
            a = self.array
            ans = reduce(np.bitwise_xor, [a[x] for x in self.neighbours(i,j)])
            if ans % 2 == 0:
                a[i, j] = 0
            else:
                a[i, j] = 1
                coords.append((i,j))
        return coords

    def likelihood(self):
        n = np.sum(self.array[self.qubit_mask()])
        N = np.sum(self.qubit_mask())
        p = self.p_val
        return p**n * (1-p)**(N-n)

    def has_logical_x_error(self):
        syndrome_sum = 0
        for i in range(0,self.L, 2):
            j = 0
            a = self.array
            nl, nr, nt, nb = a[i-1, j], a[(i+1)%self.L, j], a[i, j-1], a[i, (j+1)%self.L]
            n = nl ^ nr ^ nt ^ nb
            syndrome_sum = syndrome_sum ^ n
        return syndrome_sum % 2 == 1

    def generate_matching(self):
        coords = self.generate_x_syndrome()
        # find coords of x anyons
        # for each one Z flip the qubits required to 
        # connect it to (0,0)
        for I, J in coords:
            # Z flip first row up to I
            for i in range(1, I+1, 2): #know first qubit is at 1
                self.array[i, 0] = self.array[i, 0] ^ 1
            # Z flip Ith column up to J
            for j in range(I%2 + 1, J+1, 2):
                self.array[I, j] = self.array[I, j] ^ 1


    def generate_errors(self):
        n_qubits = np.sum(self.qubit_mask())
        errors = np.random.rand(n_qubits) < self.p_val
        self.array[self.qubit_mask()] = errors

    def generate_next(self, **kwargs):
        newstate = kwargs['newstate']
        if not newstate:
            s = self
        else:
            s = State(self.L, self.p_val)
            s.array = self.array.copy()
        # pick a random z site
        xs, ys = np.where(s.z_stabiliser_mask()==1)
        n = len(xs)
        r = np.random.random_integers(0, n-1)
        x, y = xs[r], ys[r]
        # apply the stabilizer at that point
        a = s.array
        L = np.shape(a)[0]
        a[x-1, y] = not a[x-1, y] # logical flip
        a[x, y-1] = not a[x, y-1] # logical flip
        a[(x+1)%L, y] = not a[(x+1)%L, y] # logical flip
        a[x, (y+1)%L] = not a[x, (y+1)%x] # logical flip
        return s
    
    def dump_s(self):
        return State.a_to_s(self.array)

    @staticmethod
    def a_to_s(array):
        def translate(elt, i, j):
            if (i+j)%2 == 0: # X or Z stabiliser
                return "." if elt==0 else '1'
            else: # qubit
                if elt == 0:
                    return "."
                elif elt == 1:
                    return "Z"
                elif elt == 2:
                    return "X"
                else:
                    return "Y"
        return '\n'.join([" ".join([translate(elt, i, j) for i, elt in enumerate(row)]) for j, row in enumerate(array)])

    @staticmethod
    def s_to_a(string):
        def t(elt):
            if elt == ".":
                return 0
            elif elt == '1':
                return 1
            elif elt == 'Z':
                return 1
            elif elt == 'X':
                return 2
            else : #elt == 'Y'
                return 3
        return np.array([[t(elt) for elt in s.split(" ")] for s in string.strip().split('\n')])




