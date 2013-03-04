import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import numpy as np
import full_conj_classes as fcc
import getopt

if __name__ == "__main__":
    n = int(sys.argv[1])

    if n not in range(256):

        print('Must have 0<= n <256')
        sys.exit(2)
    start = n * 256
    stop = (n + 1) * 256
    filename = './fcc_6_{}.csv'.format(n)

    fcc.write_hist_file(6, filename, start, stop)

