from ..lib import state as st
import numpy as np
import unittest
#from ..lib.conjugacy_classes import bitsum


s1 = { 'array' : """
. Z . . . Z
Z . . . Z .
. . . . . Z
. . . . . .
. Z . . . .
Z . Z . . . """, 'logical_x_error': False}

s2 = { 'array': """
. . 1 Z 1 . 1 . 1 . . .
. . . . Z . Z . . . Z .
. Z 1 . . Z 1 Z 1 . 1 Z
. . . . . . . . . . Z .
1 . . Z . . 1 Z . . . Z
. . Z . Z . . . Z . . .
1 Z . . . Z . . . Z . Z
Z . . . . . Z . . . . .
. . 1 . 1 . . Z . Z . .
Z . Z . Z . . . . . Z .
1 Z 1 Z 1 . 1 Z . . 1 Z
. . . . Z . . . Z . Z . """, 
'matching': """
. . . Z . . . Z . . . .
. . . . . . . . . . . .
. . . Z . Z . . . Z . .
. . . . . . . . . . . .
. Z . Z . Z . . . . . .
. . . . . . . . . . . .
. . . . . . . . . . . .
Z . . . . . . . . . . .
. . . Z . . . . . . . .
Z . . . . . . . . . . .
. . . Z . . . Z . Z . .
. . . . . . . . . . . .
""", 'has_logical_x_error': True}

class TestUniformToricState(unittest.TestCase):
    def assertArrayEqual(self, a, b):
        self.assertTrue(np.array_equal(a,b))

    def testQubits(self):
        s = st.ToricLattice(6)
        x = s.qubit_line('row', inbetween='z')
        self.assertArrayEqual(x, [0,0,0])
    def testMultiply(self):
        pass

    def testToStringAndBack(self):
        s = st.ToricLattice(6)
        s.change_class(13) # make state non-trivial
        s.syndrome()
        string = st.ToricLattice.a_to_s(s._array)
        a2 = st.ToricLattice.s_to_a(string)
        self.assertTrue(np.array_equal(s._array, a2))


    def testCompareApply(self):
        s = st.ToricLattice(16)
        vert_z, hor_z = st.ToricLattice.VERT_Z, st.ToricLattice.HOR_Z
        vert_x, hor_x = st.ToricLattice.VERT_X, st.ToricLattice.HOR_X
        s_orig = s.copy()
        self.assertEqual(0, st.ToricLattice.compare(s, s_orig))

        s.change_class(vert_z)
        self.assertEqual(vert_z, st.ToricLattice.compare(s, s_orig))
        s.change_class(vert_z)
        self.assertEqual(0, st.ToricLattice.compare(s, s_orig))

        s.change_class(hor_z)
        self.assertEqual(hor_z, st.ToricLattice.compare(s, s_orig))
        s.change_class(hor_z)
        self.assertEqual(0, st.ToricLattice.compare(s, s_orig))

        s.change_class(vert_x)
        self.assertEqual(vert_x, st.ToricLattice.compare(s, s_orig))
        s.change_class(vert_x)
        self.assertEqual(0, st.ToricLattice.compare(s, s_orig))

        s.change_class(vert_z + hor_x)
        self.assertEqual(vert_z + hor_x, st.ToricLattice.compare(s, s_orig))
        s.change_class(vert_z + hor_x)
        self.assertEqual(0, st.ToricLattice.compare(s, s_orig))

    def testFromSyndrome(self):

        string = """
. . 1 Z 1 . 1 . 1 . . .
. . . . Z . Z . . . Z .
. Z 1 . . Z 1 Z 1 . 1 Z
. . . . . . . . . . Z .
1 . . Z . . 1 Z . . . Z
. . Z . Z . . . Z . . .
1 Z . . . Z . . . Z . Z
Z . . . . . Z . . . . .
. . 1 . 1 . . Z . Z . .
Z . Z . Z . . . . . Z .
1 Z 1 Z 1 . 1 Z . . 1 Z
. . . . Z . . . Z . Z . """
        s = st.ToricLattice.from_string(12, string)
        synd = s.syndrome()
        s2 = st.ToricLattice.from_syndrome(12, synd)
        state_string = s2.to_s()

        string2 = """
. . . Z . . . Z . . . .
. . . . . . . . . . . .
. . . Z . Z . . . Z . .
. . . . . . . . . . . .
. Z . Z . Z . . . . . .
. . . . . . . . . . . .
. . . . . . . . . . . .
Z . . . . . . . . . . .
. . . Z . . . . . . . .
Z . . . . . . . . . . .
. . . Z . . . Z . Z . .
. . . . . . . . . . . .""".strip() # remove initial \n
        self.assertEqual(string2, state_string)

    def testNErrors(self):
        s = st.ToricLattice(16)
        # do some stuff
        s.change_class(5)
        s.apply_stabiliser(0,2)
        # check everything agrees
        n = s.n_errors()
        self.assertEqual(n, s.count_errors())
        s2 = s.copy()
        self.assertEqual(n, s2.n_errors())
        # chuck in some errors
        for i,j in s.qubit_indices:
            s.flip_qubit(i,j,flip_type = (i+j)%4)
        # check we still agree
        n = s.n_errors()
        self.assertEqual(n, s.count_errors())


