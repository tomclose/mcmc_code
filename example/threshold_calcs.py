import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import state as st
import chain as ch
import getopt

def run_trials(size=6, prob=0.1, n_steps = 5000, n_trials = 1000, verbose=True):
    successes = failures = 0

    for i in range(n_trials):
        s_orig = st.UniformToricState(size, prob)

        s_orig.generate_errors()
        synd = s_orig.syndrome()

        # assume that logical horizontal/vertical errors
        # are independent ??, so that the overall threshold
        # is just the threshold for correcting one sort
        # of error
        vert_z = s_orig.VERT_Z
        vert_x = s_orig.VERT_X
        # forget the original error confic
        s = st.UniformToricState.from_syndrome(size, prob, synd)
        conj_class = st.UniformToricState.compare(s, s_orig)
        # project onto the vertical error classes
        vert_class = conj_class & (vert_z + vert_x)

        p1 = ch.path(s.copy().change_class(vert_class)) # the right one
        p2 = ch.path(s.copy().change_class(vert_z ^ vert_class)) # with 
        p3 = ch.path(s.copy().change_class(vert_x ^ vert_class)) # with 
        p4 = ch.path(s.copy().change_class((vert_x + vert_z) ^ vert_class)) # with 

        paths = [p1, p2, p3, p4]

        ps = ch.path_set(*[ch.in_jumps_of(n_steps/2, ch.average_err(p)) for p in paths])

        # take two steps, so that the average is found
        # over the second half of the chain
        ps.next()
        res = ps.next()
        totals = [n for (n, prop, count) in res]
        if np.argmin(totals) == 0:
            successes +=1
            if verbose:
                print("success")
        else:
            failures +=1
            if verbose:
                print("fail")
    return (successes, failures)


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:n:p:L:v", ["steps=", "n_trials="])
    except getopt.GetoptError:
        usage = """
Useage:
python threshold_calcs.py -L 6 -p 0.3 --steps 1000 --n_trials 50
        """
        print(usage)
        sys.exit(2)
    # defaults:
    size, p, steps, n_trials, verbose = 6, 0.1, 1000, 20, False
    # from options:
    for opt, arg in opts:
        if opt=="-L":
            size = int(arg)
        elif opt in ("-s", "--steps"):
            steps = int(arg)
        elif opt in ("-n", "--n_trials"):
            n_trials = int(arg)
        elif opt=="-p":
            p = float(arg)
            print(p)
        elif opt == "-v":
            verbose = True
    print(run_trials(size, p, steps, n_trials, verbose=verbose))

