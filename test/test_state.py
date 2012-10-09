from ..lib.state import *
import numpy as np
import unittest

class TestState(unittest.TestCase):
    def testGenerateErrors(self):
        s = State(6, 0)
        s.generate_errors
        self.assertEqual(0, np.sum(s.array[s.qubit_mask()]))


