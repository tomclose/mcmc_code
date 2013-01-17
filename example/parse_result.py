import re
import sys
import getopt
import glob
import pickle
import numpy as np

def parse_file(file_name):
    rstart = re.compile('^Run number')
    rline = re.compile('^n: (\d+), totals: \[(\d+), (\d+), (\d+), (\d+)\], changes: \[(\d+), (\d+), (\d+), (\d+)\]')
    with open(file_name) as f:
        try:
            f.next()
            l = f.next() #get second line
        except StopIteration:
            print("File {} is not a results file".format(file_name))
            return None
        if re.match('Size: (\d+), Prob: (\d+\.\d+), n_steps: (\d+), n_trials: (\d+)', l):
            m = re.match('Size: (\d+), Prob: (\d+\.\d+), n_steps: (\d+), n_trials: (\d+)', l)
            size, p, n_steps, n_trials = m.groups()
        else:
            print("{} is not a results file".format(file_name))
            return None

        res = np.zeros((int(n_trials), 4 * int(n_steps)/1000))
        i = -1 # so the first run increases it to 0

        for line in f:
            if rstart.match(line):
                i = i+1
                j = 0
            elif rline.match(line):
                m = rline.match(line)
                n, t0, t1, t2, t3, c0, c1, c2, c3 = [int(x) for x in m.groups()]
                res[i, 4*j:4*(j+1)] = t0, t1, t2, t3
                j += 1
    return (res, p, size)


def successes(x):
    return (x[:, 0] < x[:, 1]) * (x[:, 0] < x[:, 2]) * (x[:, 0] < x[:, 3])

def results_at_time(x, n):
    return x[:, 4*n:4*(n+1)]


def gather_results(directory):
    results = []
    names = glob.glob(directory)
    for file_name in names:
        print(file_name)
        parsed = parse_file(file_name)
        if parsed is not None:
            (a, prob, size) = parsed
            prob = float(prob)
            size = int(size)
            b = successes(a[:, -4:])
            p = np.sum(b)*1.0/len(b)
            results.append([prob, size, p])
    return results


def load_results(file_name):
    with open(file_name) as f:
        res = pickle.load(f)
    return res

def plot_results(res):
    pass
        

if __name__ == "__main__":
    try:
        outfile = sys.argv[1]
    except getopt.GetoptError:
        usage = """parse_results.py outfile"""
    if outfile is None:
        outfile = './results.pickle'
    res = gather_results("./*")
    with open(outfile, 'w') as f:
        pickle.dump(res, f)




"""
Run number: 100
Errors: [0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0
 0 3 3 0 0 3 0 0 1 0 0 2 0]
Syndrome: [(0, 0), (1, 5), (3, 5), (5, 5), (6, 6), (6, 8), (7, 5), (7, 9), (8, 0), (8, 4), (8, 8), (9, 7)]
Class: 1
n: 1000, totals: [9588, 17136, 12267, 13506], changes: [71, 175, 116, 123]
n: 2000, totals: [19025, 28414, 23530, 27575], changes: [122, 257, 206, 242]
n: 3000, totals: [28394, 38476, 35677, 40255], changes: [171, 317, 302, 321]
n: 4000, totals: [38595, 48272, 46632, 53821], changes: [233, 382, 399, 426]
n: 5000, totals: [46948, 58609, 57728, 66035], changes: [283, 445, 480, 497]
n: 6000, totals: [57160, 70068, 68402, 79157], changes: [360, 515, 544, 589]
n: 7000, totals: [65721, 80898, 80880, 94371], changes: [408, 599, 638, 714]
n: 8000, totals: [74983, 92556, 96571, 109937], changes: [460, 675, 796, 833]
n: 9000, totals: [83797, 103919, 107489, 124190], changes: [520, 737, 878, 93
n: 10000, totals: [91922, 114160, 118105, 140651], changes: [565, 790, 969, 1044]
Result: 1
"""
