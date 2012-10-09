import sys
import os
# this path stuff is so retarded
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import state

s = state.State(6, 0.3)
s.generate_errors()
s.generate_x_syndrome()
s.show()
for i in range(10):
    s.generate_next()
    s.generate_x_syndrome()
    s.show()
