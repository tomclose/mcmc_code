from ..lib.state import *
import numpy as np
import unittest
from ..lib.conjugacy_classes import bitsum


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
    def testGenerateErrors(self):
        s = UniformToricState(6, 0)
        s.generate_errors()
        self.assertEqual(0, np.sum(s.array[zip(*s.qubit_indices)]))
    def testToStringAndBack(self):
        s = UniformToricState(6, 0.5)
        s.generate_errors()
        s.generate_syndrome()
        string = UniformToricState.a_to_s(s.array)
        a2 = UniformToricState.s_to_a(string)
        self.assertTrue(np.array_equal(s.array, a2))
    def testXError(self):
        s = UniformToricState(6, 0)
        s.array = UniformToricState.s_to_a(s1['array'])
        self.assertFalse(s.measure_hor_x_loop())
    def testS2(self):
        s = UniformToricState(12, 0.3)
        s.array = UniformToricState.s_to_a(s2['array'])

        #check x syndrome generation works
        s.generate_syndrome()
        self.assertArrayEqual(s.array, UniformToricState.s_to_a(s2['array']))

        # check matching works
        s.generate_matching()
        self.assertArrayEqual(s.matching, UniformToricState.s_to_a(s2['matching']))

        # check that has_logical_error works
        self.assertTrue(s.measure_hor_x_loop)
    def testGenerateNext(self):
        s = UniformToricState(16, 0.4)
        s.generate_errors()
        x_synd = s.generate_syndrome()
        matching = s.generate_matching()
        for i in range(5):
            s.generate_next()
            # syndrome shouldn't change
            self.assertArrayEqual(x_synd, s.generate_syndrome())
        # when apply matching it should lie in codespace
        s.array = s.array ^ matching
        self.assertArrayEqual([], s.generate_syndrome())

    def testToN(self):
        """ ToN should only work on arrays of 'small' dimesnion (<6?)"""
        s = UniformToricState(6, 0.4)
        # 0 errors should always map to 0
        self.assertEqual(s.to_n(), 0)
        # find stabiliser representation
        stab_n = s.apply_stabiliser(0,0).to_n()
        # check stabiliser rep looks sensible
        self.assertEqual(bitsum(stab_n), 4)
        # do commutativity consistency check
        s.generate_errors()
        n1 = s.to_n()^stab_n # convert then apply stab
        n2 = s.apply_stabiliser(0,0).to_n() # apply stab then convert
        self.assertEqual(n1, n2)

    def testMeasureApply(self):
        s = ToricLattice(16)
        self.assertFalse(s.measure_hor_x_loop())
        self.assertFalse(s.measure_vert_x_loop())
        self.assertFalse(s.measure_hor_z_loop())
        self.assertFalse(s.measure_vert_z_loop())

        s.add_hor_x_loop()
        self.assertTrue(s.measure_vert_x_loop())
        s.add_hor_x_loop()
        self.assertFalse(s.measure_vert_x_loop())

        s.add_hor_z_loop()
        self.assertTrue(s.measure_vert_z_loop())
        s.add_hor_z_loop()
        self.assertFalse(s.measure_vert_z_loop())

        s.add_vert_x_loop()
        self.assertTrue(s.measure_hor_x_loop())
        s.add_vert_x_loop()
        self.assertFalse(s.measure_hor_x_loop())

        s.add_vert_z_loop()
        self.assertTrue(s.measure_hor_z_loop())
        s.add_vert_z_loop()
        self.assertFalse(s.measure_hor_z_loop())
