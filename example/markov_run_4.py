import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import state as st
import chain as ch

def run(size=6, prob=0.1):
    s = st.UniformToricState(size, prob)
    s.generate_errors()
    s.generate_matching()
    # make no error the right answer
    hor_z = s.measure_hor_z_loop()
    vert_z = s.measure_vert_z_loop()

    # forget the original error config
    s.array[:] = s.matching[:]
    # and make the 1st conjugacy class the correct answer
    if vert_z:
        print("Applying hor z")
        s.add_hor_z_loop()
    if hor_z:
        print("Applying vert z")
        s.add_vert_z_loop()


    p = ch.path(s)
    pv = ch.path(s.copy().add_hor_z_loop())
    ph = ch.path(s.copy().add_vert_z_loop())
    phv = ch.path(s.copy().add_vert_z_loop().add_hor_z_loop())

    paths = [p, pv, ph , phv]


    ps = ch.path_set(*[ch.n_jumps(ch.average_err(p), 1000) for p in paths])

    old = np.array([0, 0, 0, 0, 0, 0, 0, 0])
    for a in ps:
        new = np.array([n for pair in a for n in pair])
        res = (new-old)/1000.0
        print('(%.2f, %.2f) (%.2f, %.2f) (%.2f, %.2f) (%.2f, %.2f)'% tuple(res))
        old = new

