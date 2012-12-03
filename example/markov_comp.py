import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import state as st
import chain as ch
import conjugacy_classes as cc

prob = 0.1
size = 6

ccs = cc.ConjugacyClasses(st.ToricLattice(size))

def go():
    s = st.UniformToricState(6, prob)
    s.generate_errors()
    s.generate_matching()
    # make no error the right answer
    if s.measure_hor_x():
        s.apply_vert_z()
    if s.measure_vert_x():
        s.apply_hor_z()
    if s.measure_hor_x() or s.measure_vert_x():
        raise "oops something went wrong with the logical ops"

    # forget the original error config
    s.array[:] = s.matching[:]


    i, v, h, vh = ccs.lookup_hist(s.to_n())

    av_i = cc.ConjugacyClasses.hist_row_av_n(i, prob)
    av_v = cc.ConjugacyClasses.hist_row_av_n(v, prob)
    av_h = cc.ConjugacyClasses.hist_row_av_n(h, prob)
    av_vh = cc.ConjugacyClasses.hist_row_av_n(vh, prob)

    p = ch.path(s)
    pv = ch.path(s.copy().apply_vert_z())
    ph = ch.path(s.copy().apply_hor_z())
    phv = ch.path(s.copy().apply_vert_z().apply_hor_z())

    paths = [p, pv, ph , phv]

    def average_err(p):
        count = 0
        total_n = 0
        for current, next_s, changed in p:
            count+=1
            total_n += current.n_errors()
            yield total_n*1.0/count

    ps = ch.path_set(*[ch.n_jumps(average_err(p), 1000) for p in paths])

    for a in ps:
        print(av_i, av_v, av_h, av_vh)
        print(a)
