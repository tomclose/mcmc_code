from ..lib import state as st
import numpy as np
import unittest
from ..lib import conjugacy_classes as cc

class TestToN(unittest.TestCase):
    def assertArrayEqual(self, a, b):
        self.assertTrue(np.array_equal(a,b))

    def testToN(self):
        """ ToN should only work on arrays of 'small' dimesnion (<6?)"""
        s = st.UniformToricState(6, 0.4)
        # 0 errors should always map to 0
        self.assertEqual(cc.to_n(s), 0)
        # find stabiliser representation
        stab_n = cc.to_n(s.apply_stabiliser(0,0))
        # check stabiliser rep looks sensible
        self.assertEqual(cc.bitsum(stab_n), 4)
        # do commutativity consistency check
        s.generate_just_z_errors()
        n1 = cc.to_n(s)^stab_n # convert then apply stab
        n2 = cc.to_n(s.apply_stabiliser(0,0)) # apply stab then convert
        self.assertEqual(n1, n2)
