import numpy as np
#import state

class Chain:
    def __init__(self, s0):
        self.states = [s0]
        self.errors = [s0.has_logical_x_error()]

    def advance(self, n):
        for i in range(n):
            self.step()

    def current_state(self):
        return self.states[-1]

    def step(self):
        current = self.current_state()
        new = current.generate_next(newstate=True)
        p = current.likelihood()
        q = new.likelihood()
        # r = q/p not robust, instead of X < r do p*X < q
        if p * np.random.rand() < q: # includes the r>1 case
            self.states.append(new)
            self.errors.append(new.has_logical_x_error())
            return True
        else:
            self.states.append(current)
            self.errors.append(current.has_logical_x_error())
            return False

    def state(self,i):
        return self.states[i]



