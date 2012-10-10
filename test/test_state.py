from ..lib.state import *
import numpy as np
import unittest


s1 = { 'array' : """
1 Z 1 . 1 Z
Z . . . Z .
1 . . . 1 Z
. . . . . .
1 Z 1 . . .
Z . Z . . . """, 'no': 5}


class TestState(unittest.TestCase):
    def testGenerateErrors(self):
        s = State(6, 0)
        s.generate_errors()
        self.assertEqual(0, np.sum(s.array[s.qubit_mask()]))
    def testToStringAndBack(self):
        s = State(6, 0.5)
        s.generate_errors()
        s.generate_x_syndrome()
        string = State.a_to_s(s.array)
        a2 = State.s_to_a(string)
        self.assertTrue(np.array_equal(s.array, a2))
        

