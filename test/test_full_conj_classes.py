from ..lib import state as st
import numpy as np
import unittest
from ..lib import full_conj_classes as fcc

def to_int(string):
    ans = 0
    l = 0
    for x in string:
        x = int(x)
        if x not in [0,1,2,3]:
            raise RuntimeError('to_int passed {}. Should be 0, 1, 2, or 3'.format(x))
        ans += (x << 2*l)
        l += 1
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

    def testBareStabiliserActions(self):
        # X 0 X 0
        # 0 Z 0 Z
        # X 0 X 0
        # 0 Z 0 Z
        actions = [ to_int('11100010'),
                    to_int('11010001'),
                    to_int('20222000'),
                    to_int('02220200'),
                    to_int('00101110'),
                    to_int('00011101'),
                    to_int('20002022'),
                    to_int('02000222') ]
        self.assertSameElts(actions, fcc.bare_stabiliser_actions(4))
