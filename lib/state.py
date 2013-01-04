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
    VERT_Z = 1
    HOR_Z = 2
    VERT_X = 4
    HOR_X = 8
    Z = 1
    X = 2
    Y = 3
    def __init__(self, L):
        self.L = L
        self._array = np.zeros((L,L), dtype='uint8')
        self.x_i_indices = range(0, L, 2)
        self.x_j_indices = range(0, L, 2)
        self.z_i_indices = range(1, L, 2)
        self.z_j_indices = range(1, L, 2)
        self.x_stabiliser_indices = [(i,j) for i in self.x_i_indices for j in self.x_j_indices]
        self.z_stabiliser_indices = [(i,j) for i in self.z_i_indices for j in self.z_j_indices]
        self.qubit_indices = [(i,j) for j in range(0, self.L) for i in range((j+1)%2, self.L, 2)]
        self.stabiliser_indices = [(i,j) for j in range(0, self.L) for i in range(j%2, self.L, 2)]
        self._n_errors = 0
        self._syndrome = None

    def flip_qubit(self, i, j, flip_type):
        val = self._array[i,j]
        new_val = val ^ flip_type
        if val == 0:
            if new_val > 0:
                self._n_errors += 1
        else: # val > 0
            if new_val == 0:
                self._n_errors -= 1
        self._array[i,j] = new_val

    def qubit(self, i, j):
        return self._array[i,j]

    def qubit_line(self, direction, inbetween):
        return self._array[zip(*self.qubit_line_ij(direction, inbetween))]

    def qubit_line_ij(self, direction, inbetween):
        """ Usage:
                qubit_line_ij('row', inbetween='x')
        """
        if direction == 'row':
            if inbetween == 'z':
                qs = [(i, self.z_j_indices[0]) for i in self.x_i_indices]
            elif inbetween == 'x':
                qs = [(i, self.x_j_indices[0]) for i in self.z_i_indices]
            else:
                raise RuntimeError("Inbetween must be 'x' or 'z'", inbetween)
        elif direction == 'col':
            if inbetween == 'z':
                qs = [(self.z_i_indices[0], j) for j in self.x_j_indices]
            elif inbetween == 'x':
                qs = [(self.x_i_indices[0], j) for j in self.z_j_indices]
            else:
                raise RuntimeError("Inbetween must be 'x' or 'z'", inbetween)
        else:
            raise RuntimeError("Invalid direction (should be 'row' or 'col')", direction)
        return qs


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
        return self._n_errors

    # count errors is relatively expensive, so only call when necessary
    def count_errors(self):
        self._n_errors = int(np.sum(self._array[zip(*self.qubit_indices)]>0))
        return self._n_errors




    # Actions
    # =======
    def syndrome(s, refresh=False):
        if s._syndrome is None or refresh:
            coords = []
            for i,j in s.stabiliser_indices:
                if s.measure_stabiliser(i,j):
                    s._array[i,j] = 1
                    coords.append((i,j))
                else:
                    s._array[i,j] = 0
            s._syndrome = coords
        return s._syndrome

    def apply_stabiliser(s, x, y):
        site_type = s.site_type(x, y)
        # site_type = 1 for an X stabiliser site,
        #             2 for a Z
        #             0 for qubit
        # NB passing in a qubit site will silently have no effect
        for i,j in s.neighbours(x,y):
            s.flip_qubit(i,j, site_type)
        return s

    def measure_stabiliser(s, x, y):
        return bool(s.site_type(x,y) & reduce(np.bitwise_xor, [s._array[c] for c in s.neighbours(x,y)]))


    def copy(self):
        s = self.__class__(self.L)
        self.copy_onto(s)
        return s

    def copy_onto(self, other_state):
        other_state._array = self._array.copy()
        other_state._n_errors = self.n_errors()
        return other_state

    def change_class(s, change):
        if change & s.VERT_Z:
            for i, j in s.qubit_line_ij('row', inbetween='x'):
                s.flip_qubit(i, j, s.Z)
        if change & s.HOR_Z:
            for i, j in s.qubit_line_ij('col', inbetween='x'):
                s.flip_qubit(i, j, s.Z)
        if change & s.VERT_X:
            for i, j in s.qubit_line_ij('row', inbetween='z'):
                s.flip_qubit(i, j, s.X)
        if change & s.HOR_X:
            for i, j in s.qubit_line_ij('col', inbetween='z'):
                s.flip_qubit(i, j, s.X)
        return s

    @staticmethod
    def multiply(s1, s2):
        s = s1.copy()
        for i,j in s.qubit_indices:
            s.flip_qubit(i, j, s2.qubit(i,j))
        return s

    @classmethod
    def compare(cls, s1, s2):
        """ Returns the number of the syndrome 
            xHor zHor xVert zVert
            i.e 4 = zHor
                5 = zHor, zVert
                15 = zHor, zVert, xHor, xVert
        """
        s = cls.multiply(s1, s2) # is self the class here
        result = 0
        result += cls.VERT_Z if reduce(lambda v, q: q^v, s.qubit_line('col', inbetween='z')) & 1 else 0
        result += cls.HOR_Z  if reduce(lambda v, q: q^v, s.qubit_line('row', inbetween='z')) & 1 else 0
        result += cls.VERT_X if reduce(lambda v, q: q^v, s.qubit_line('col', inbetween='x')) & 2 else 0
        result += cls.HOR_X  if reduce(lambda v, q: q^v, s.qubit_line('row', inbetween='x')) & 2 else 0
        return result
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

    @classmethod
    def from_syndrome(cls, L, syndrome, new_obj=None):
        if new_obj is not None:
            s = new_obj # to allow descendants to pass in 
                        # a blank instance of themselves
        else:
            s = cls(L)
        for i, j in syndrome:
            flip_type = s.site_type(i, j) # 1 for X, 2 for Z
            while(j > 1): # move to the left
                # jump a qubit an flip it
                s.flip_qubit(i, j-1, flip_type)
                j -= 2
            while(i > 1): # move to the top
                # jump a qubit and flip it
                s.flip_qubit(i-1, j, flip_type)
                i -= 2
        s.count_errors()
        return s

    # Displaying
    # ==========
    def to_s(self):
        return self.__class__.a_to_s(self._array)


    def show(self):
        fig = plt.gcf()# or plt.figure()
        fig.clf()
        ax = fig.add_subplot(111)
        # Decid how to show array:
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
        #for x, y in np.argwhere(self.matching > 0):
            #circ = plt.Circle((y, x), radius = 0.3)
            # it seems that matplotlib transposes the coords
            # of an imshow such that the coords where 
            # a[i,j] are shown are actually (j, i)
            # this makes sense, as for arrays j corresponds
            # to horizontal movement, whereas the convention
            # for graph coords is (x=horiz, y=vert)
            #ax.add_patch(circ)
        plt.show() # in case not already drawn
        plt.draw() # refresh if drawn
        return show_array

    @classmethod
    def from_string(cls, n, string):
        s = cls(n)
        s._array = cls.s_to_a(string)
        s.count_errors()
        return s

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
        return x**diff

    def generate_errors(self):
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
        self.count_errors()

    def copy(self):
        s = self.__class__(self.L, self.p)
        self.copy_onto(s)
        return s

    def generate_next(s):
        # pick a random z site
        n = len(s.stabiliser_indices)
        r = np.random.random_integers(0, n-1)
        x, y = s.stabiliser_indices[r]
        # apply the stabilizer at that point
        s.apply_stabiliser(x,y)

    @classmethod
    def from_syndrome(cls, L, p, syndrome):
        s = cls(L, p)
        return ToricLattice.from_syndrome(L, syndrome, s)


class ZUniformToricState(UniformToricState):
    def generate_errors(self):
        n_qubits = len(self.qubit_indices)
        errors = np.random.rand(n_qubits) < self.p
        self._array[zip(*self.qubit_indices)] = errors
        self.count_errors()

    def generate_next(s):
        # pick a random z site
        n = len(s.x_stabiliser_indices)
        r = np.random.random_integers(0, n-1)
        x, y = s.x_stabiliser_indices[r]
        # apply the stabilizer at that point
        s.apply_stabiliser(x,y)

class XUniformToricState(UniformToricState):
    def generate_errors(self):
        n_qubits = len(self.qubit_indices)
        errors = np.random.rand(n_qubits) < self.p
        self._array[zip(*self.qubit_indices)] = 2*errors
        self.count_errors()

    def generate_next(s):
        # pick a random x site
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
