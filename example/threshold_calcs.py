import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import state as st
import chain as ch

def run_trials(size=6, prob=0.1, n_steps = 5000, n_trials = 1000):
    successes = failures = 0
    s = st.UniformToricState(size, prob)

    for i in range(n_trials):
        s.generate_errors()
        m = s.generate_matching()
        hor_z = s.measure_hor_z_loop()
        vert_z = s.measure_vert_z_loop()

        # forget the original error confic
        s.array[:] = m[:]

        p = ch.path(s)
        pv = ch.path(s.copy().add_hor_z_loop())
        ph = ch.path(s.copy().add_vert_z_loop())
        phv = ch.path(s.copy().add_vert_z_loop().add_hor_z_loop())

        if not hor_z and not vert_z:
            # identity is the right answer
            paths = [p, pv, ph, phv]
        elif hor_z and not vert_z:
            # we have a vert_z_loop
            paths = [pv, p, ph, phv]
        elif not hor_z and vert_z:
            # we have a hor_z_loop
            paths = [ph, p, pv, phv]
        else:
            paths = [phv, p, pv, ph]

        ps = ch.path_set(*[ch.n_jumps(ch.average_err(p), n_steps/2) for p in paths])
        
        # take two steps, so that the average is found
        # over the second half of the chain
        ps.next()
        res = ps.next()
        totals = [n for (n, prop) in res]
        if np.argmin(totals) == 0:
            successes +=1
            print("success")
        else:
            failures +=1
            print("fail")
    return (successes, failures)


        

