import numpy as np
import state as st
import chain as ch
import matplotlib.pyplot as plt

s = st.UniformToricState(6, 0.3)
p = ch.path(s)

x = np.array([x for x in ch.aggregate(ch.first_n(p, 100))])

plt.plot(x)
