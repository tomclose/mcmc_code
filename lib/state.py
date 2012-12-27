import numpy as np
import matplotlib.pyplot as plt


class ToricLattice:
    """ State is reponsible for holding the configuration
    of a toric lattice - both errors and stabiliser states.
    
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

    Measure_vert_z detects any logical vertical z errors, by measuring for zs along a vertical z row.
    Apply_vert_z applies a vertical z error, by flipping along a horizontal x row
    
    X 1 X 1 X 1     <-- this is a vertical Z error
    . Z . Z . Z         it can't be seen by the Xs
    X . X . X .
    . Z . Z . Z
    X . X . X .
    . Z . Z . Z
    
    """
    def __init__(self, L):
        self.L = L
        self._array = np.zeros((L,L), dtype='uint8')
        self._matching = None
        self.x_i_indices = range(0, L, 2)
        self.x_j_indices = range(0, L, 2)
        self.z_i_indices = range(1, L, 2)
        self.z_j_indices = range(1, L, 2)
        self.x_stabiliser_indices = [(i,j) for i in self.x_i_indices for j in self.x_j_indices]
        self.z_stabiliser_indices = [(i,j) for i in self.z_i_indices for j in self.z_j_indices]
        self.qubit_indices = [(i,j) for j in range(0, self.L) for i in range((j+1)%2, self.L, 2)]
        self._n_errors = None

    def flip_qubit(self, i, j, flip_type='Z'):
        val = self._array[i,j]
        if flip_type == 'Z':
            new_val = val ^ 1
        elif flip_type == 'X':
            new_val = val ^ 2
        elif flip_type == 'Y':
            new_val = val ^ 3
        else:
            raise RuntimeError("Invalid flip_type: ", flip_type)
        self._array[i,j] = new_val


    @property
    def matching(self):
        if self._matching is None:
            self.reset_matching()
        return self._matching

    def neighbours(self, i, j):
        return [(i-1, j), (i, j-1), ((i+1)%self.L, j), (i, (j+1)%self.L)]
    def site_type(self, i, j):
        """ Returns:
                0 - for a qubit
                1 - for an X plaquette (as they detect z errors)
                2 - for a Z plaquette
        """
        i_mod2, j_mod2 = i%2, j%2
        if i_mod2==0 and j_mod2==0:
            return 1
        elif i_mod2==1 and j_mod2==1:
            return 2
        else:
            return 0

    # Querying
    # ========
    
    def n_errors(self):
        # cast as an int, otherwise it returns a uint8, which
        # leads to all types of problems if you try to do 
        # arithmetic (e.g. 113-115 = 383498534198239348)
        #if self._n_errors is None:
        self._n_errors = int(np.sum(self._array[zip(*self.qubit_indices)]))
        return self._n_errors

    """ 
    Measure_vert_z detects any logical vertical z errors, by measuring for zs along a vertical z row.
    Apply_vert_z applies a vertical z error, by flipping along a horizontal x row
    
    X 1 X 1 X 1     <-- this is a vertical Z error
    . Z . Z . Z         it can't be seen by the Xs
    X . X . X .
    . Z . Z . Z
    X . X . X .
    . Z . Z . Z
    
    """
    def has_hor_z(s, m = None):
        m = m or s.matching
        error_sum = 0
        # Z-measure the 1s along a Z-col
        j = s.z_j_indices[0]
        for i in s.x_i_indices:
            n = s._array[i, j]^m[i,j]
            error_sum = error_sum ^ n
        return error_sum & 1 == 1
    def has_hor_x(s, m = None):
        m = m or s.matching
        error_sum = 0
        # X-measure the 2s along a X-col
        j = s.x_j_indices[0]
        for i in s.z_i_indices:
            n = s._array[i, j]^m[i,j]
            error_sum = error_sum ^ n
        return error_sum & 2 == 2
    def has_vert_x(s, m = None):
        m = m or s.matching
        error_sum = 0
        # X-measure the 2s along a X-row
        i = s.x_i_indices[0]
        for j in s.z_j_indices:
            n = s._array[i, j]^m[i,j]
            error_sum = error_sum ^ n
        return error_sum & 2 == 2
    def has_vert_z(s, m=None):
        m = m or s.matching
        error_sum = 0
        # Z-measure the 1st along a Z-row
        i = s.z_i_indices[0]
        for j in s.x_j_indices:
            n = s._array[i, j]^m[i,j]
            error_sum = error_sum ^ n
        return error_sum & 1 == 1

    def apply_hor_x(s):
        s._n_errors = None # reset error_count
        # X-flip the qubits on a Z-row
        i = s.z_i_indices[0]
        for j in s.x_j_indices:
            s.flip_qubit(i, j, 'X')
        return s
    def apply_vert_x(s):
        s._n_errors = None # reset error_count
        # X-flip the qubits on a Z-column
        j = s.z_j_indices[0]
        for i in s.x_i_indices:
            s.flip_qubit(i, j, 'X')
        return s
    def apply_hor_z(s):
        s._n_errors = None # reset error_count
        # Z-flip the qubits on an X-row
        i = s.x_i_indices[0]
        for j in s.z_j_indices:
            s.flip_qubit(i, j, 'Z')
        return s
    def apply_vert_z(s):
        s._n_errors = None # reset error_count
        # Z-flip the qubits on an X-column
        j = s.x_j_indices[0]
        for i in s.z_i_indices:
            s.flip_qubit(i, j, 'Z')
        return s
    # Actions
    # =======
    def generate_syndrome(s):
        s.reset_syndrome()
        xs = s.generate_x_syndrome()
        zs = s.generate_z_syndrome()
        return xs + zs

    def reset_syndrome(s):
        s._array[zip(*s.z_stabiliser_indices)] = 0
        s._array[zip(*s.x_stabiliser_indices)] = 0

    def generate_z_syndrome(s):
        coords = []
        for i, j in s.x_stabiliser_indices:
            if s.measure_stabiliser(i,j):
                s._array[i, j] = 1
                coords.append((i,j))
            else:
                s._array[i,j] = 0
        return coords

    def generate_x_syndrome(s):
        coords = []
        for i, j in s.z_stabiliser_indices:
            if s.measure_stabiliser(i,j):
                s._array[i, j] = 1
                coords.append((i,j))
            else:
                s._array[i,j] = 0
        return coords

    def generate_matching(self):
        self.reset_matching()
        self.generate_z_matching()
        self.generate_x_matching()
        return self._matching
        
    def reset_matching(self):
        self._n_errors = None
        self._matching = np.zeros(np.shape(self._array), dtype='uint8')

    def generate_z_matching(self):
        m = self._matching
        self.generate_z_syndrome() # do twice so don't change state *yuk*
        coords = self.generate_z_syndrome()
        for I, J in coords:
            # Z flip first row up to I
            for i in range(1, I+1, 2): #know first qubit is at 1
                m[i, 0] = m[i, 0] ^ 1
            # Z flip Ith column up to J
            for j in range((I+1)%2, J+1, 2):
                m[I, j] = m[I, j] ^ 1
        return m

    def generate_x_matching(self):
        m = self._matching
        self.generate_x_syndrome() # do twice to avoid state change ** yuk **
        coords = self.generate_x_syndrome()
        for I, J in coords: #I, J will be odd
            # X flip second row up to I
            for i in range(2, I, 2): #know first qubit is at 1
                m[i, 1] = m[i, 1] ^ 2
            # Z flip Ith column up to J
            for j in range(2, J, 2):
                m[I, j] = m[I, j] ^ 2
        return m
    def apply_stabiliser(s, x, y):
        site_type = s.site_type(x, y)

        if site_type == 1:   # X plaquette
            flip_type = 'Z'
        elif site_type == 2: # Z plaquette
            flip_type = 'X'
        else:                # a qubit
            raise RuntimeError("Not a stabiliser site", x, y)

        for i,j in s.neighbours(x,y):
            s.flip_qubit(i,j, flip_type)

        return s

    def measure_stabiliser(s, x, y):
        return bool(s.site_type(x,y) & reduce(np.bitwise_xor, [s._array[c] for c in s.neighbours(x,y)]))


    def copy(self):
        s = self.__class__(self.L)
        self.copy_onto(s)
        return s

    def copy_onto(self, other_state):
        other_state._array = self._array.copy()
        # note that matching is NOT copied - it's a reference
        # this is fine - each time we call generate_matching()
        #                a new one is created
        other_state.matching = self.matching
        return other_state


    # Displaying
    # ==========
    def dump_s(self):
        return __class__.a_to_s(self._array)


    def show(self):
        fig = plt.gcf()# or plt.figure()
        fig.clf()
        ax = fig.add_subplot(111)
        # Decide how to show array:
        #  - 0 is a non-firing X stabiliser
        #  - 1 is a non-firing Z stabiliser
        #  - 2 is a non-error qubit
        #  - 6 is an error qubit
        #  - 4 is a firing stabiliser
        show_array = np.zeros((self.L, self.L), dtype='int')
        a = self._array
        for (i,j) in self.x_stabiliser_indices:
            show_array[i,j] = 4 if a[i,j]==1 else 0
        for (i,j) in self.z_stabiliser_indices:
            show_array[i,j] = 4 if a[i,j]==1 else 1
        for (i,j) in self.qubit_indices:
            if a[i,j] ==1:
                show_array[i,j] = 6
            elif a[i,j] ==2:
                show_array[i,j] = 7
            elif a[i,j] ==3:
                show_array[i,j] = 8
            else:
                show_array[i,j] = 2
        cax = ax.imshow(show_array, interpolation='nearest', vmax=8)
        cbar = plt.colorbar(cax, ticks=[0, 1, 2, 4, 6, 7, 8])
        cbar.ax.set_yticklabels(['X stab', 'Z stab','qubit', 'firing X stab', 'Z error qubit',  'X error qubit',  'Y error qubit'])
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


class UniformToricState(ToricLattice):
    def __init__(self, L, p):
        ToricLattice.__init__(self, L)
        self.p = p

    def likelihood(self):
        n = self.n_errors()
        #N = len(self.qubit_indices)
        p = self.p/(1-self.p)
        return p**n #* (1-p)**(N-n)

    def relative_prob(self, s2):
        # returns p(s2)/p(self)
        n = self.n_errors()
        n2 = s2.n_errors()
        x = self.p/(1-self.p)
        diff = n2 - n
        #print(n2, n, self.p,x, diff, x**diff )

        return x**diff

    def generate_just_z_errors(self):
        self._n_errors = None # reset error_count
        n_qubits = len(self.qubit_indices)
        errors = np.random.rand(n_qubits) < self.p
        self._array[zip(*self.qubit_indices)] = errors

    def generate_just_x_errors(self):
        self._n_errors = None # reset error_count
        n_qubits = len(self.qubit_indices)
        errors = np.random.rand(n_qubits) < self.p
        self._array[zip(*self.qubit_indices)] = 2*errors

    def generate_errors(self):
        self._n_errors = None # reset error_count
        n_qubits = len(self.qubit_indices)
        r = np.random.rand(n_qubits)
        def map_to_error(x):
            if x < self.p/3:
                return 3
            elif x < 2*self.p/3:
                return 2
            elif x < self.p:
                return 1
            else:
                return 0
        self._array[zip(*self.qubit_indices)] = [map_to_error(x) for x in r]

    def copy(self):
        s = self.__class__(self.L, self.p)
        self.copy_onto(s)
        return s

    def generate_next(s):
        s._n_errors = None # reset error_count
        # pick a random z site
        n = len(s.z_stabiliser_indices)
        r = np.random.random_integers(0, n-1)
        x, y = s.z_stabiliser_indices[r]
        # apply the stabilizer at that point
        s.apply_stabiliser(x,y)






class HotState(ToricLattice):
    def generate_next(s):
        if np.random.rand() < 0.5: 
            # make a logical x error,
            # by flipping a z row
            for j in range(1, s.L, 2):
                i = 0
                s._array[i, j] = s._array[i,j] ^ 1
        # pick a random z site
        for x, y in s.z_stabiliser_indices:
            if np.random.rand() < 0.5:
                # apply the stabilizer at that point
                s.apply_stabiliser(x,y)
