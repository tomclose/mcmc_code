import sys
import os
# this path stuff is so retarded
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import state
import chain

s = state.State(6, 0.5)
s.generate_errors()
s.generate_matching()
s.generate_x_syndrome()
s.show()

c = chain.Chain(s)
for i in range(20):
    print(c.step())
    c.current_state().show()

