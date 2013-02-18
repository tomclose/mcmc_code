from ..lib import state as st
import numpy as np
import unittest
from ..lib import y_conjugacy_classes as ycc

class TestYConjClasses(unittest.TestCase):
    def assertArrayEqual(self, a, b):
        self.assertTrue(np.array_equal(a,b))
    def assertSameElts(self, a, b):
        a.sort()
        b.sort()
        self.assertArrayEqual(a, b)

    def testBitSum(self):
        self.assertEqual(ycc.bitsum(5), 2)
        self.assertEqual(ycc.bitsum(32), 1)
        self.assertEqual(ycc.bitsum(31), 5)
        self.assertEqual(ycc.bitsum(7), 3)

    def testToNAndSetState(self):
        s = st.ToricLattice(4)
        self.assertEqual(ycc.errors_to_n(s), 0)
        ycc.set_state_errors(s, 13)
        self.assertEqual(ycc.errors_to_n(s), 13)
        ycc.set_state_errors(s, 7)
        self.assertEqual(ycc.errors_to_n(s), 7)

    def testSyndromeToN(self):
        # x 0 * 0
        # 0 * 3 *
        # x 0 * 0   syndrome => 00101110
        # 0 z 0 z   error => 00001000
        s = st.ToricLattice(4)
        self.assertEqual(ycc.syndrome_to_n(s), 0)
        ycc.set_state_errors(s, int('00001000', 2))
        s.syndrome(refresh=True)
        self.assertEqual(ycc.syndrome_to_n(s), int('00101110', 2))

    def testClassToN(self):
        # x 0 x 3
        # 0 z 3 z   class => zHor + zVert + xHor + xVert
        # x 3 x 0   syndrome => 00000000
        # 3 z 0 z   error => 01011010
        s = st.ToricLattice(4)
        self.assertEqual(ycc.n_to_class(s, 0), 0)
        self.assertEqual(ycc.n_to_class(s, int('01011010', 2)), s.VERT_Z + s.HOR_Z + s.VERT_X + s.HOR_X)

    def testNToSyndrome(self):
        # x 0 * 0
        # 0 * 3 *
        # x 0 * 0   syndrome => 00101110
        # 0 z 0 z   error => 00001000
        s = st.ToricLattice(4)
        self.assertEqual(ycc.n_to_syndrome(s, 0), 0)
        self.assertEqual(ycc.n_to_syndrome(s, int('00001000', 2)), int('00101110', 2))
