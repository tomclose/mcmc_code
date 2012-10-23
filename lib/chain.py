import numpy as np
#import state

class Chain:
    def __init__(self, s0):
        self.current_state = s0
        self.next_state = s0.copy()
        self.error_types = s0.error_types

    def advance(self, n):
        for i in range(n):
            self.step()

    def step(self, next_state = None):
        current = self.current_state
        if next_state is not None:
            new = next_state.copy_onto(self.next_state)
        else: # copy the current state and generate a new one
            new = current.copy_onto(self.next_state)
            new.generate_next()
        p = current.likelihood()
        q = new.likelihood()
        # r = q/p not robust, instead of X < r do p*X < q
        if p * np.random.rand() < q: # includes the r>1 case
            self.next_state.copy_onto(self.current_state)
            return True
        else:
            return False


