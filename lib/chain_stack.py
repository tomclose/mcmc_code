import numpy as np
import state
import chain

def create_chainstack(s0, n_chains, n_steps):
    p = s0.p_val
    probabilities = [p + (0.5-p)/(n_chains-1)*i for i in range(n_chains)]
    states = [s0.copy() for i in range(n_chains-1)]
    hot_state = state.HotState(s0.L, 0.5)
    s0.copy_onto(hot_state)
    states.append(hot_state)
    for p, s in zip(probabilities, states):
        s.p_val = p
    # set errors on first state
    # probably need a mechanism to set errors on all states from first
    chains = [chain.Chain(s) for s in states]
    return ChainStack(chains, n_steps)

class ChainStack:
    def __init__(self, chains, steps):
        self.chains = chains
        # put the original chain at 0
        self.step_counter = np.zeros(len(chains), dtype='int')
        self.n_chains = len(chains)
        self.n_chain_steps = steps 
        self.error_types = chains[0].error_types

    def twist(self):
        """Tries a swap operation between the chains.

        Returns a list of the indices of any chains that took on a new state
        """
        cs = self.chains
        successes = []
        # if 5 chains (0-4) I will try to copy 4 onto 3, ... , 1 onto 0
        # Note: range(3, 0, -1) = [3,2,1]
        for i in range(self.n_chains-1, 0, -1):
            if self.metropolis_says_yes(cs[i].current_state, cs[i-1].current_state):
                cs[i].current_state.copy_onto(cs[i-1].current_state)
                successes.append(i-1)
        return successes

    def step(self):
        n_steps = self.n_chain_steps
        twist_successes = []

        # if other chains aren't fully advanced
        # run each of the other chains forward n steps
        for i in range(1, self.n_chains):
            n_to_go = n_steps - self.step_counter[i]
            if n_to_go > 0:
                self.chains[i].advance(n_to_go)
            self.step_counter[i] = self.n_chain_steps

        # if main chain isn't fully advanced
        if self.step_counter[0] < n_steps:
            # run the principal chain forward one step
            self.chains[0].step()
            self.step_counter[0] += 1

        # if we've now just done the final step
        if self.step_counter[0] == n_steps:
            # do the twist
            twist_successes = self.twist()
            # reset counter
            self.step_counter *= 0

        s = self.chains[0].current_state
        step_results = {'error_count': s.n_errors(), 'error_type':s.logical_error()}

        return  {'state_details': step_results, 'twists': twist_successes}


    @staticmethod
    def metropolis_says_yes(state1, state2):
        p = state1.likelihood()
        q = state2.likelihood()
        r = q/p
        return np.random.rand() < r # includes the r>1 case
