from ..lib import state as st
import numpy as np
import unittest
from ..lib import full_conj_classes as fcc

def to_int(string):
    ans = 0
    for x in string:
        x = int(x)
        if x not in [0,1,2,3]:
            raise RuntimeError('to_int passed {}. Should be 0, 1, 2, or 3'.format(x))
        ans = (ans << 2) + x
    return ans

class TestYConjClasses(unittest.TestCase):
    def assertArrayEqual(self, a, b):
        self.assertTrue(np.array_equal(a,b))
    def assertSameElts(self, a, b):
        a.sort()
        b.sort()
        self.assertArrayEqual(a, b)

    def testToNAndSetState(self):
        s = st.ToricLattice(4)
        self.assertEqual(fcc.errors_to_n(s), 0)
        fcc.set_state_errors(s, 13)
        self.assertEqual(fcc.errors_to_n(s), 13)
        fcc.set_state_errors(s, 7)
        self.assertEqual(fcc.errors_to_n(s), 7)

    def testStabActionForSite(self):
        # X 0 X 0
        # 0 Z 0 Z
        # X 0 X 0  Zs breed 1s, Xs breed 2s
        # 0 Z 0 Z
        s = st.ToricLattice(4)
        results = [(0, to_int('02000222'))]
        for site, action in results:
            self.assertEqual(action, fcc.stab_action_for_site(s, site))

    def testBareStabiliserActions(self):
        # X 0 X 0
        # 0 Z 0 Z
        # X 0 X 0  Zs breed 1s, Xs breed 2s
        # 0 Z 0 Z
        actions = [ to_int('02000222'),
                    to_int('20002022'),
                    to_int('00011101'),
                    to_int('00101110'),
                    to_int('02220200'),
                    to_int('11010001')]
        self.assertSameElts(actions, fcc.bare_stabiliser_actions(4))

    def testCompositeStabiliserAction(self):
        bsa = fcc.bare_stabiliser_actions(4)
        self.assertEqual(0, fcc.composite_stabiliser_action(4, 0, bsa))
        # * 0 X 0
        # 0 Z 0 * error = to_int('02101332')
        # X 0 X 0 stabiliser = int('001010', 2) [read ^-- for Zs starting at 
        # 0 Z 0 Z                                 second last, then same for X]
        self.assertEqual(to_int('02101332'),
                         fcc.composite_stabiliser_action(4,
                                                         int('010001', 2),
                                                         bsa))
    def testCompositeStabiliserActions(self):
        self.assertSameElts(fcc.composite_stabiliser_actions(2), [0])
        # sanity check for n = 4
        acts4 = fcc.composite_stabiliser_actions(4)
        self.assertEqual(len(acts4), 2**6)

    def testClassCandidates(self):
        # X 0
        # 0 Z
        self.assertSameElts(fcc.class_candidates(2), range(16))


    def testErrorCandidateForSyndrome(self):
        state = st.ToricLattice(2)
        self.assertEqual(0, fcc.error_candidate_for_syndrome(0, 2))
        s = st.ToricLattice(4)


    def testSyndromeCandidates(self):
        self.assertSameElts([0], [ x for x in fcc.syndrome_candidates(2)])

    def testErrorCount(self):
        self.assertEqual(fcc.error_count(int('032103', 4)), 4)

    def testSyndromeClassOrbits(self):
        x2 = np.array([o for s,c,o in fcc.syndrome_class_orbits(2)])
        x2 = x2.reshape((-1,))
        x2.sort()
        self.assertArrayEqual(x2, range(16))
        x4 = np.array([o for s,c,o in fcc.syndrome_class_orbits(4)])
        x4 = x4.reshape((-1,))
        x4.sort()
        self.assertArrayEqual(x4, range(4**8))

    def testErrorDist(self):
        orb = [ 0, 1, 2, 3, 4, 5, 6, 7]
        self.assertEqual(fcc.error_dist(orb), {0:1, 1:4, 2:3})

    def testHistRow(self):
        orb = range(8)
        self.assertArrayEqual(fcc.hist_row(orb, 3), [1, 4, 3, 0])

    def testHistArray(self):
        ha2 = fcc.hist_array(2)[:, 2:]
        count2 = np.sum(ha2, axis=0)
        # coeffs in (1+3x)**2
        self.assertArrayEqual([1, 6, 9], count2)
        ha4 = fcc.hist_array(4)[:, 2:]
        count4 = np.sum(ha4, axis=0)
        # coeffs in (1+3x)**8
        self.assertArrayEqual([1, 24, 252, 1512, 5670, 13608, 20412, 17496, 6561], count4)

    def testClassProbabilities(self):
        ha2 = fcc.hist_array(2)
        ha4 = fcc.hist_array(4)
        self.assertAlmostEqual(1, sum(fcc.class_probabilities(ha2, 0.01)[:,2]))
        self.assertAlmostEqual(1, sum(fcc.class_probabilities(ha2, 0.22)[:,2]))
        self.assertAlmostEqual(1, sum(fcc.class_probabilities(ha2, 0.50)[:,2]))
        self.assertAlmostEqual(1, sum(fcc.class_probabilities(ha2, 0.88)[:,2]))
        self.assertAlmostEqual(1, sum(fcc.class_probabilities(ha4, 0.01)[:,2]))
        self.assertAlmostEqual(1, sum(fcc.class_probabilities(ha4, 0.22)[:,2]))
        self.assertAlmostEqual(1, sum(fcc.class_probabilities(ha4, 0.50)[:,2]))
        self.assertAlmostEqual(1, sum(fcc.class_probabilities(ha4, 0.88)[:,2]))

    def testSuccessProbability(self):
        pass

    def testBinaryCombinations(self):
        self.assertSameElts(fcc.binary_combinations(1, 0), [0])
        self.assertSameElts(fcc.binary_combinations(1, 1), [1])
        self.assertSameElts(fcc.binary_combinations(2, 1), [1, 2])
        self.assertSameElts(fcc.binary_combinations(2, 2), [3])
        self.assertSameElts(fcc.binary_combinations(4, 2), [3, 5, 9, 6, 10, 12])

    def testXSyndsNearZero(self):
        # (0)000 => (1)000, (0)001, (0)010, (0)100
        # (1)000 => (0)000, (1)001, (1)010, (1)100
        self.assertSameElts(fcc.x_synds_near_zero(3, 1), [0, 1, 2, 4])


    def testSyndsNearZero(self):
        # (0)00(0)00, 2 => (1)01(0)00, (1)10(0)00, (0)11(0)00, (0)00(1)01, (0)00(1)10, (0)00(0)11
        self.assertSameElts(fcc.synds_near_zero(2, 2, 0, 0), [4, 8, 12, 1, 2, 3])
