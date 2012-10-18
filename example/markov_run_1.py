import chain_stack
import state
import markov_run

s = state.State(6, 0.3)
s.generate_errors()
s.generate_matching()
cs = chain_stack.create_chainstack(s, 3, 10)
m = markov_run.MarkovRun(cs)

for i in range(10000):
    m.step()
