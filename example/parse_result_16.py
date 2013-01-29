import re
import sys
import getopt
import glob
import pickle

def parse_file(f):
    try:
        f.next()
        l = f.next() #get second line
    except StopIteration:
        print("File {} is not a results file".format(f.name))
        raise StopIteration
    if re.match('Size: (\d+), Prob: (\d+\.\d+), n_steps: (\d+), n_trials: (\d+)', l):
        m = re.match('Size: (\d+), Prob: (\d+\.\d+), n_steps: (\d+), n_trials: (\d+)', l)
        size, p, n_steps, n_trials = m.groups()
        p = float(p)
        size = int(size)
        yield (size, p, contentParser(f))
    else:
        print("{} is not a results file".format(f.name))
        raise StopIteration

class contentParser:
    """
    Usage:
        for expt in contentParser(f):
            for n, vals in expt:
                print("After {n} iterations values were {vals}")
    """
    def __init__(self, file_handle):
        self.file_it = iter(file_handle)
        self._need_new_line = True
        self.rstart = re.compile('^Run number')
        self.rline = re.compile('^n: (\d+), totals: \[(.*?)\].*changes: \[(.*?)\]')
    def __iter__(self):
        return self
    def next_line(self):
        if self._need_new_line:
            self._current_line = next(self.file_it)
            self._need_new_line = False
        return self._current_line
    def line_used(self):
        # mark current line as stale
        # gets around pulling a new line before necessary
        # and the associated risk of hitting the
        # StopIteration exception too soon
        self._need_new_line = True
    def next(self):
        while(1):
            if self.rstart.match(self.next_line()):
                self.line_used()
                return self._lines()
            else:
                self.line_used()
    def _lines(self):
        line_match = self.rline.search(self.next_line())
        # first move to the next matching line
        while line_match is None:
            self.line_used()
            line_match = self.rline.search(self.next_line())
        # then yield each of the block of matching lines
        while line_match:
            n = int(line_match.groups()[0])
            res = line_match.groups()[1]
            res = [int(r) for r in res.split(",")]
            yield [n, res]
            self.line_used()
            line_match = self.rline.search(self.next_line())

def successes(x):
    return (x[:, 0] < x[:, 1]) * (x[:, 0] < x[:, 2]) * (x[:, 0] < x[:, 3]) * (x[:, 0] < x[:, 4]) * (x[:, 0] < x[:, 5]) * (x[:, 0] < x[:, 6]) * (x[:, 0] < x[:, 7]) * (x[:, 0] < x[:, 8]) * (x[:, 0] < x[:, 9]) * (x[:, 0] < x[:, 10]) * (x[:, 0] < x[:, 11]) * (x[:, 0] < x[:, 12]) * (x[:, 0] < x[:, 13]) * (x[:, 0] < x[:, 14]) * (x[:, 0] < x[:, 15])

def results_at_time(x, n):
    return x[:, 16*n:16*(n+1)]

def gather_results(pattern):
    results = []
    names = glob.glob(pattern)
    for file_name in names:
        print(file_name)
        with open(file_name) as f:
            for (size, prob, runs) in parse_file(f):
                n = 0
                n_success = 0
                for run in runs:
                    lines = [l for l in run]
                    if len(lines) ==0:
                        print("No lines in file {}".format(file_name))
                        break
                    step_n, final_totals = lines[-1]
                    f = final_totals[0]
                    success = reduce(lambda x, y: x and y,
                                        [f <= v for v in final_totals[1:]],
                                        True)
                    n += 1
                    if success: n_success += 1
                if n==0:
                    p = 0
                else:
                    p = n_success*1.0/n
                results.append([prob, size, n_success, n])
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
