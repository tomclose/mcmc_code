import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import state as st
import chain as ch
import conjugacy_classes as cc

prob = 0.1

s = st.UniformToricState(6, prob)
ccs = cc.ConjugacyClasses(s)
s.generate_errors()

p = ch.path(s)


i, v, h, vh = ccs.lookup_hist(s.to_n())

av_i = cc.ConjugacyClasses.hist_row_av_n(i, prob)

total_n = 0
count = 0


for current, next_s, changed in p:
    count+=1
    total_n += current.n_errors()
    if count%100000 == 0:
        print(total_n*1.0/count, av_i)
        raw_input()



