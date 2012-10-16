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
        self.matching = np.zeros((L,L), dtype='bool')
        self.p_val = p
        self.next_state = None

    def show(self):
        fig = plt.gcf()# or plt.figure()
        fig.clf()
        ax = fig.add_subplot(111)
        # Decide how to show array:
        #  - 0 is a non-firing X stabiliser
        #  - 1 is a non-firing Z stabiliser
        #  - 2 is a non-error qubit
        #  - 5 is an error qubit
        #  - 6 is a firing stabiliser
        show_array = np.zeros((self.L, self.L), dtype='int')
        show_array[self.z_stabiliser_mask()] = 1 # Z stab = 1
        show_array[self.qubit_mask()] = 2
        show_array[np.where(self.x_stabiliser_mask()*self.array == 1)] = 4
        show_array[np.where(self.qubit_mask()*self.array == 1)] = 6
        cax = ax.imshow(show_array, interpolation='nearest', vmax=6)
        cbar = plt.colorbar(cax, ticks=[0, 1, 2, 4, 6])
        cbar.ax.set_yticklabels(['X stab', 'Z stab','qubit', 'firing X stab', 'error qubit'])
        for x, y in np.argwhere(self.matching > 0):
            circ = plt.Circle((y, x), radius = 0.3)
            # it seems that matplotlib transposes the coords
            # of an imshow such that the coords where 
            # a[i,j] are shown are actually (j, i)
            # this makes sense, as for arrays j corresponds
            # to horizontal movement, whereas the convention
            # for graph coords is (x=horiz, y=vert)
            ax.add_patch(circ)
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
        error_sum = 0
        # measure out the qubits along a Z column
        for i in range(0, self.L, 2):
            j = 1
            n = self.array[i, j]^self.matching[i,j]
            error_sum = error_sum ^ n
        return error_sum & 1 == 1

    def generate_matching(self):
        coords = self.generate_x_syndrome()
        self.matching = np.zeros(np.shape(self.array), dtype='bool')
        # find coords of x anyons
        # for each one Z flip the qubits required to 
        # connect it to (0,0)
        m = self.matching
        for I, J in coords:
            # Z flip first row up to I
            for i in range(1, I+1, 2): #know first qubit is at 1
                m[i, 0] = m[i, 0] ^ 1
            # Z flip Ith column up to J
            for j in range((I+1)%2, J+1, 2):
                m[I, j] = m[I, j] ^ 1
        return m


    def generate_errors(self):
        n_qubits = np.sum(self.qubit_mask())
        errors = np.random.rand(n_qubits) < self.p_val
        self.array[self.qubit_mask()] = errors
    
    def set_next_state(self, state):
        self.next_state = state

    def generate_next(self, **kwargs):
        newstate = kwargs.get('newstate')
        if self.next_state:
            s = self.next_state
            if not newstate:
                self.array = s.array.copy()
                s = self
            self.next_state = None
        else:
            if not newstate:
                s = self
            else:
                s = State(self.L, self.p_val)
                s.array = self.array.copy()
                s.matching = self.matching
            # pick a random z site
            n = len(s.z_stabiliser_indices())
            r = np.random.random_integers(0, n-1)
            x, y = s.z_stabiliser_indices()[r]
            # apply the stabilizer at that point
            for i,j in s.neighbours(x,y):
                s.array[i,j] = s.array[i,j] ^ 1
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




