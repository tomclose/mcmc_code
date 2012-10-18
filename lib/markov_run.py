import numpy as np

class MarkovRun:
    def __init__(self, chain_stack):
        self.chain_stack = chain_stack
        self.reset()
    
    def reset(self):
        self.twist_count = np.zeros(self.chain_stack.n_chains-1)
        self.results = np.zeros((100,6))
        self.capacity = 100
        self.size = 0
        self.t = 0 # should basically always be the same as size
        self.total_error_count = 0 # keeps track of the running total of errors
                                   # in the 0 state at each stage
        self.error_types = self.chain_stack.error_types
        self.error_type_count = np.zeros(len(self.error_types))

    def append_results(self, row):
        if self.size == self.capacity:
            # expand array size
            self.capacity *= 4
            new_results = np.zeros((self.capacity, 6))
            new_results[:self.size, 6] = self.results[:,6]
            self.results = new_results

        self.results[self.size, 6] = row[:]
        self.size += 1

    def step(self):
        res = self.chain_stack.step()

        n_chains = self.chain_stack.n_chains
        twists = res['twists']
        tcount = self.twist_count
        for n in sorted(twists,reverse=True):
            # want to increase position n only if n-1 is equal
            # and n+1 is ahead:
            # 3,3,4 => 3,4,4
            # 2,3,4 => 2,3,4
            # 3,3,3 => 3,3,3

            # say we have 4 chains
            # then n could be 0, 1, 2
            behind_check = n==n_chains-2 or tcount[n+1]-tcount[n] ==1
            ahead_check = n==0 or tcount[n-1] == tcount[n]
            if behind_check and ahead_check:
                tcount[n] += 1

        # results array:
        #  [tops0, total_error_count, errors_types.....]
        tops0 = self.twist_count[0]
        self.total_error_count += res['state_details']['error_count']
        logical_error_index = self.error_types.index(res['state_details']['error_type'])
        self.error_type_count[logical_error_index] += 1

        print(tops0, self.total_error_count, self.error_type_count, self.twist_count)
        
