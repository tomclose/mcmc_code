from ..lib import state as st
import numpy as np
import unittest
from ..lib import conjugacy_classes as cc

class TestToN(unittest.TestCase):
    def assertArrayEqual(self, a, b):
        self.assertTrue(np.array_equal(a,b))
    def assertSameElts(self, a, b):
        a.sort()
        b.sort()
        self.assertArrayEqual(a, b)

    def testToN(self):
        """ ToN should only work on arrays of 'small' dimesnion (<6?)"""
        s = st.ZUniformToricState(6, 0.4)
        # 0 errors should always map to 0
        self.assertEqual(cc.to_n(s), 0)
        # find stabiliser representation
        stab_n = cc.to_n(s.apply_stabiliser(1,1))
        # check stabiliser rep looks sensible
        self.assertEqual(cc.bitsum(stab_n), 4)
        # do commutativity consistency check
        s.generate_errors()
        n1 = cc.to_n(s)^stab_n # convert then apply stab
        n2 = cc.to_n(s.apply_stabiliser(1,1)) # apply stab then convert
        self.assertEqual(n1, n2)

    def testSetState(self):
        s = st.ToricLattice(4)
        # set up a random state
        for i, j in s.qubit_indices:
            if np.random.rand() > 0.5:
                s.flip_qubit(i, j, s.Z)
        old_qubits = [s.qubit(i,j) for i, j in s.qubit_indices]
        n = cc.to_n(s)
        # mess up state
        for i, j in s.qubit_indices:
            s.flip_qubit(i, j, s.Z)
        # reset state
        cc.set_state(s, n)
        new_qubits = [s.qubit(i,j) for i, j in s.qubit_indices]
        self.assertEqual(old_qubits, new_qubits)

    def testStabGenEffect(self):
        # X 1 X 0
        # 1 Z 1 Z
        # X 1 X 0
        # 0 Z 0 Z # => 00011101 = 
        state = st.ToricLattice(4)
        x = cc.stab_gen_effect(state, 1, 1)
        self.assertEqual(x, int('00011101', 2))

    def testStabGens(self):
        state = st.ToricLattice(4)
        stab_gens = cc.stab_gens(state)
        actual_stab_gens = [int('00011101', 2),
                            int('00101110', 2),
                            int('11010001', 2),
                            int('11100010', 2)]
        self.assertSameElts(stab_gens, actual_stab_gens)

    def testStabEffect(self):
        state = st.ToricLattice(4)
        stab_gens = cc.stab_gens(state)
        # this is the 1001 = 9th stab
        # it will have the same effect as 0110
        # * 1 X 1
        # 1 Z 1 Z
        # X 1 * 1
        # 1 Z 1 Z # => 11111111 = 
        x = cc.stab_effect(stab_gens, int('1001',2))
        self.assertEqual(x, int('11111111', 2))

    def testStabilisers(self):
        state = st.ToricLattice(4)
        # there are 8 independent stabilisers
        # note that stabilisers are labelled starting
        # at top left, as are qubits
        # stabiliser 1:
        # X 1 X 0
        # 1 * 1 Z
        # X 1 X 0
        # 0 Z 0 Z # => 00011101
        actual_stabs = [ 0,
                        int('00011101', 2), #1
                        int('00101110', 2), #2
                        int('00110011', 2), #3
                        int('11010001', 2), #4
                        int('11001100', 2), #5
                        int('11111111', 2), #6
                        int('11100010', 2)] #7
        stabs = cc.stabilisers(state)
        self.assertSameElts(stabs, actual_stabs)

    def testOrbit(self):
        state = st.ToricLattice(4)
        stabs = cc.stabilisers(state)
        # the orbit of 0 should be the stabilisers
        self.assertSameElts(cc.orbit(stabs, 0), stabs)
        # now look at the orbit of 1
        actual_orb = [int('00000001', 2), # apply 0
                      int('00011100', 2), # apply 1
                      int('00101111', 2), # apply 2
                      int('00110010', 2), # apply 3
                      int('11010000', 2), # apply 4
                      int('11001101', 2), # apply 5
                      int('11111110', 2), # apply 6
                      int('11100011', 2)] # apply 7
        self.assertSameElts(cc.orbit(stabs, 1), actual_orb)

    def testConjClassGens(self):
        # X 0 X 0
        # 1 Z 0 Z
        # X 0 X 0
        # 1 Z 0 Z # => 01000100
        state = st.ToricLattice(4)
        actual_ccs = [ 0,
                       int('01000100', 2), # hor Z
                       int('00000011', 2), # vert z
                       int('01000111', 2)] # both
        self.assertSameElts(cc.conj_class_gens(state), actual_ccs)

    def testMatching(self):
        # X 0 X 0
        # 1 Z 0 Z
        # X 1 * 0
        # 0 Z 0 Z # => synd = 1000, state = 00010100
        state = st.ToricLattice(4)
        self.assertEqual(int('00010100', 2), cc.matching(16, state))

    def testErrorDist(self):
        state = st.ToricLattice(4)
        stabs = cc.stabilisers(state)
        # test the 0 syndrome
        o1 = cc.orbit(stabs, 0)
        d1 = {0:1, 4:6, 8:1}
        self.assertEqual(cc.error_dist(o1), d1)
        # test the 1 syndrome
        o2 = cc.orbit(stabs, 1)
        d2 = {1:1, 3:3, 5:3, 7:1}
        self.assertEqual(cc.error_dist(o2), d2)

    def testHistRow(self):
        error_dist = {0:4, 3:2, 5:1}
        correct_answer = [4, 0, 0, 2, 0, 1]
        answer = cc.hist_row(error_dist, 5)
        self.assertEqual(correct_answer, answer)

    def testOverlapArray(self):
        a2 = cc.overlap_array(2)
        r2 = [[0], [1]]
        self.assertArrayEqual(a2, r2)
        r4 = [[0, 2, 2, 2, 2, 2, 2, 4],
              [2, 0, 2, 2, 2, 2, 4, 2],
              [2, 2, 0, 2, 2, 4, 2, 2],
              [2, 2, 2, 0, 4, 2, 2, 2],
              [2, 2, 2, 4, 0, 2, 2, 2],
              [2, 2, 4, 2, 2, 0, 2, 2],
              [2, 4, 2, 2, 2, 2, 0, 2],
              [4, 2, 2, 2, 2, 2, 2, 0],
              [1, 1, 1, 3, 1, 3, 3, 3],
              [1, 1, 3, 1, 3, 1, 3, 3],
              [1, 3, 1, 1, 3, 3, 1, 3],
              [3, 1, 1, 1, 3, 3, 3, 1],
              [1, 3, 3, 3, 1, 1, 1, 3],
              [3, 1, 3, 3, 1, 1, 3, 1],
              [3, 3, 1, 3, 1, 3, 1, 1],
              [3, 3, 3, 1, 3, 1, 1, 1]]
        a4 = cc.overlap_array(4)
        self.assertArrayEqual(a4, r4)

    def testSyndIndexedBy(self):
        self.assertEqual(int('1001', 2), cc.synd_indexed_by(int('0001', 2), 4))
        self.assertEqual(int('0011', 2), cc.synd_indexed_by(int('0011', 2), 4))
        self.assertEqual(int('1011', 2), cc.synd_indexed_by(int('1011', 2), 4))
        self.assertEqual(int('0001', 2), cc.synd_indexed_by(int('1001', 2), 4))

    def testSyndromeProbs(self):
        n = 4
        p = 0.15
        pp = cc.syndrome_probs(n, p)
        n_stabs = (n/2)**2
        print(n_stabs)
        norm_odd = (1 - (1 - 2*p)**n_stabs)/2
        norm_even = (1 + (1 - 2*p)**n_stabs)/2
        n_rows = 2**n_stabs
        results = np.ones(n_rows)
        results[0:n_rows/2] *= norm_even
        results[n_rows/2:] *= norm_odd
        totals = np.sum(pp, axis=1)
        np.testing.assert_allclose(totals, results)


    def testNoisyProb(self):
        pass
