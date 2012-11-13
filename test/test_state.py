from ..lib.state import *
import numpy as np
import unittest


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

class TestState(unittest.TestCase):
    def assertArrayEqual(self, a, b):
        self.assertTrue(np.array_equal(a,b))
    def testGenerateErrors(self):
        s = State(6, 0)
        s.generate_errors()
        self.assertEqual(0, np.sum(s.array[zip(*s.qubit_indices())]))
    def testToStringAndBack(self):
        s = State(6, 0.5)
        s.generate_errors()
        s.generate_x_syndrome()
        string = State.a_to_s(s.array)
        a2 = State.s_to_a(string)
        self.assertTrue(np.array_equal(s.array, a2))
    def testXError(self):
        s = State(6, 0.5)
        s.array = State.s_to_a(s1['array'])
        self.assertFalse(s.has_logical_x_error())
    def testS2(self):
        s = State(12, 0.3)
        s.array = State.s_to_a(s2['array'])

        #check x syndrome generation works
        s.generate_x_syndrome()
        self.assertArrayEqual(s.array, State.s_to_a(s2['array']))

        # check matching works
        s.generate_matching()
        self.assertArrayEqual(s.matching, State.s_to_a(s2['matching']))

        # check that has_logical_error works
        self.assertTrue(s.has_logical_x_error)
    def testGenerateNext(self):
        s = State(16, 0.4)
        s.generate_errors()
        x_synd = s.generate_x_syndrome()
        matching = s.generate_matching()
        for i in range(5):
            s.generate_next()
            # syndrome shouldn't change
            self.assertArrayEqual(x_synd, s.generate_x_syndrome())
        # when apply matching it should lie in codespace
        s.array = s.array ^ matching
        self.assertArrayEqual([], s.generate_x_syndrome())

    def testToN(self):
        s = State(16, 0.4)
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


