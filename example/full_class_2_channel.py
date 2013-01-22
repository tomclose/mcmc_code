import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import state as st
import chain as ch
import getopt

""" Note that size here is not L as in state.py, L = 2 * size """

def run_trials(size=6, prob=0.1, n_steps = 5000, n_trials = 1000, verbose=True):
    successes = failures = 0
    if verbose: print("Size: {}, Prob: {}, n_steps: {}, n_trials: {}".format(size,
        prob, n_steps, n_trials))

    for i in range(n_trials):
        s_orig = st.UniformToricState(2 * size, prob)

        if verbose: print("Run number: {}".format(i+1))
        s_orig.generate_errors()
        if verbose: print("Errors: {}".format(s_orig.qubits()))
        synd = s_orig.syndrome()
        if verbose: print("Syndrome: {}".format(synd))

        s = st.UniformToricState.from_syndrome(2 *size, prob, synd)

        conj_class = st.UniformToricState.compare(s, s_orig)
        if verbose: print("Class: {}".format(conj_class))

        # 16 conjugacy classes, reorder so the right one is first
        conj_classes = [conj_class ^ i for i in range(16)]

        paths = [ch.path(s.copy().change_class(c)) for c in conj_classes]
        accum_paths = [ch.average_err(p) for p in paths ]
        chunked_paths = [ch.in_jumps_of(1000, p) for p in accum_paths]
        ps = ch.path_set(*chunked_paths)

        # take two steps, so that the average is found
        # over the second half of the chain
        for res in ps:
            totals = [n for (n, prop, count) in res]
            changes = [prop for (n, prop, count) in res]
            n = res[0][2]
            if verbose: print("n: {}, totals: {}, changes: {}".format(n, totals, changes))
            if n >= n_steps:
                break

        if np.argmin(totals) == 0:
            successes +=1
            if verbose: print("Result: 1")
        else:
            failures +=1
            if verbose: print("Result: 0")

        if verbose: print("\n")
    if verbose: print("Size: {}, Prob: {}, n_steps: {}, n_trials: {}".format(size,
        prob, n_steps, n_trials))
    if verbose: print("Total successes: {}".format(successes))
    if verbose: print("Total failures: {}".format(failures))
    if verbose: print("Total proportion: {}".format(successes*1.0/(successes + failures)))
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

