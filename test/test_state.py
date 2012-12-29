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

    def testQubits(self):
        s = ToricLattice(6)
        x = s.qubit_line('row', inbetween='z')
        self.assertArrayEqual(x, [0,0,0])
    def testMultiply(self):
        pass

    def testGenerateErrors(self):
        s = UniformToricState(6, 0)
        s.generate_errors()
        self.assertEqual(0, np.sum(s._array[zip(*s.qubit_indices)]))
    def testToStringAndBack(self):
        s = UniformToricState(6, 0.5)
        s.generate_errors()
        s.syndrome()
        string = UniformToricState.a_to_s(s._array)
        a2 = UniformToricState.s_to_a(string)
        self.assertTrue(np.array_equal(s._array, a2))
    def testXError(self):
        s = UniformToricState(6, 0)
        s.array = UniformToricState.s_to_a(s1['array'])
        self.assertFalse(s.has_hor_x())
    def testS2(self):
        s = UniformToricState(12, 0.3)
        s.array = UniformToricState.s_to_a(s2['array'])
        #check x syndrome generation works
        s.syndrome()
        print(s.array)
        print(UniformToricState.s_to_a(s2['array']))
        self.assertArrayEqual(s._array, UniformToricState.s_to_a(s2['array']))

        # check matching works
        s.generate_matching()
        self.assertArrayEqual(s.matching, UniformToricState.s_to_a(s2['matching']))

        # check that has_logical_error works
        self.assertTrue(s.has_hor_x)

    def testSyndromeMatchingX(self):
        s = UniformToricState(12, 0.3)
        s.generate_just_x_errors()
        s.reset_matching()
        m = s.generate_x_matching()
        s._array = s._array ^ m
        self.assertArrayEqual([], s.syndrome())

    def testSyndromeMatchingZ(self):
        s = UniformToricState(12, 0.3)
        s.reset_matching()
        s.generate_just_z_errors()
        m = s.generate_z_matching()
        s._array = s._array ^ m
        self.assertArrayEqual([], s.syndrome())

    def testGenerateNext(self):
        s = UniformToricState(16, 0.4)
        s.generate_errors()
        x_synd = s.syndrome()
        matching = s.generate_matching()
        for i in range(5):
            s.generate_next()
            # syndrome shouldn't change
            self.assertArrayEqual(x_synd, s.syndrome())
        # when apply matching it should lie in codespace
        print(s._array)
        s._array = s._array ^ matching
        print(s._array)
        self.assertArrayEqual([], s.syndrome())


    def testCompareApply(self):
        s = ToricLattice(16)
        s_orig = s.copy()
        self.assertEqual(0, ToricLattice.compare(s, s_orig))

        s.change_class(8)
        self.assertEqual(8, ToricLattice.compare(s, s_orig))
        s.change_class(8)
        self.assertEqual(0, ToricLattice.compare(s, s_orig))

        s.change_class(4)
        self.assertEqual(4, ToricLattice.compare(s, s_orig))
        s.change_class(4)
        self.assertEqual(0, ToricLattice.compare(s, s_orig))

        s.change_class(2)
        self.assertEqual(2, ToricLattice.compare(s, s_orig))
        s.change_class(2)
        self.assertEqual(0, ToricLattice.compare(s, s_orig))

        s.change_class(1)
        self.assertEqual(1, ToricLattice.compare(s, s_orig))
        s.change_class(1)
        self.assertEqual(0, ToricLattice.compare(s, s_orig))
