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
print("Logical x error: ", s.has_logical_x_error())

c = chain.Chain(s)
for i in range(20):
    c.step()
    print(c.errors[-1])
    s = c.current_state()
    s.generate_x_syndrome()
    c.current_state().show()

