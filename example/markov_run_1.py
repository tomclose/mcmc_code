import chain_stack
import state
import markov_run

s = state.State(20, 0.5)
s.generate_errors()
s.generate_matching()
cs = chain_stack.create_chainstack(s, 3, 9)
m = markov_run.MarkovRun(cs)

for i in range(1000):
    m.step()
    print(m.current_state())

