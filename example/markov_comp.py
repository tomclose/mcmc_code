import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import state as st
import chain as ch
import conjugacy_classes as cc

size = 6

ccs = cc.ConjugacyClasses(st.ToricLattice(size))

def go(prob):
    s = st.UniformToricState(6, prob)

    #s.generate_just_z_errors()
    m = s.generate_matching()
    #hor_z = s.has_hor_z()
    #vert_z = s.has_vert_z()

    # forget the original error confic
    s.array[:] = m[:]

    sv = s.copy().apply_hor_z()
    sh = s.copy().apply_vert_z()
    shv = s.copy().apply_vert_z().apply_hor_z()

    p = ch.path(s)
    pv = ch.path(sv)
    ph = ch.path(sh)
    phv = ch.path(shv)

    #if not hor_z and not vert_z:
        ## identity is the right answer
    paths = [p, pv, ph, phv]
    #elif hor_z and not vert_z:
        ## we have a vert_z_loop
        #paths = [pv, p, ph, phv]
    #elif not hor_z and vert_z:
        ## we have a hor_z_loop
        #paths = [ph, p, pv, phv]
    #else:
        #paths = [phv, p, pv, ph]


    i, v, h, vh = ccs.lookup_hist(cc.to_n(s))

    av_i = cc.ConjugacyClasses.hist_row_av_n(i, prob)
    av_v = cc.ConjugacyClasses.hist_row_av_n(v, prob)
    av_h = cc.ConjugacyClasses.hist_row_av_n(h, prob)
    av_vh = cc.ConjugacyClasses.hist_row_av_n(vh, prob)

    n_steps = 1000
    ps = ch.path_set(*[ch.n_jumps(ch.average_err(p), n_steps) for p in paths])

    for a in ps:
        print(av_i, av_v, av_h, av_vh)
        print([b[0]*1.0/b[2] for b in a])
