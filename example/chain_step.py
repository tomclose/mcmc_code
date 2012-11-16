import sys
import os
# this path stuff is so retarded
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import state
import chain

s = state.State(12, 0.3)
s.generate_errors()
s.generate_matching()
s.generate_x_syndrome()
s.show()
print("Logical x error: ", s.has_logical_x_error())

c = chain.Chain(s)
for (s, new, changed) in (c.next() for i  in range(30)):
    s.show()
    print(changed)

