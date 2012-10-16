import numpy as np
import state
import chain

def create_chainstack(L, p, n, steps):
    states = [state.State(L, p + (0.5-p)/(n-1)*i) for i in range(n)]
    # set errors on first state
    # probably need a mechanism to set errors on all states from first
    chains = [chain.Chain(s) for s in states]
    return ChainStack(chains, steps)

class ChainStack:
    def __init__(self, chains, steps):
        self.chains = chains
        # put the original chain at 0
        self.chain_switches = np.zeros(len(chains)-1)
        self.n_chains = len(chains)
        self.n_chain_steps =steps 

    def twist(self):
        cs = self.chains
        switch = np.zeros(self.n_chains-1)
        for i in range(self.n_chains-1, 0, -1):
            if self.metropolis_says_yes(cs[i].current_state(), cs[i-1].current_state()):
                cs[i-1].current_state().set_next_state(cs[i].current_state())
                switch[i-1] = 1
        self.chain_switches = np.vstack([self.chain_switches, switch])

    def step(self):
        for c in self.chains:
            c.advance(self.n_chain_steps)
        self.twist()

    def advance(self, n):
        for i in range(n):
            self.step()

    @staticmethod
    def metropolis_says_yes(state1, state2):
        p = state1.likelihood()
        q = state2.likelihood()
        r = q/p
        return np.random.rand() < r # includes the r>1 case
