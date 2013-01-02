import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import state as st
import chain as ch
import conjugacy_classes as cc

size = 6

ccs = cc.ConjugacyClasses(st.ToricLattice(size))

def go(prob):
    hor_z = st.ToricLattice.HOR_Z
    vert_z = st.ToricLattice.VERT_Z

    s_orig = st.UniformToricState(6, prob)

    s_orig.generate_just_z_errors()
    synd = s_orig.syndrome()

    s = st.UniformToricState.from_syndrome(6, prob, synd)

    sv = s.copy().change_class(vert_z)
    sh = s.copy().change_class(hor_z)
    shv = s.copy().change_class(hor_z + vert_z)

    p = ch.path(s)
    pv = ch.path(sv)
    ph = ch.path(sh)
    phv = ch.path(shv)

    paths = [p, pv, ph, phv]

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
