#!/usr/bin/env python
'''Take docking commands file as argument and check that outputs are reasonable.'''

import sys, re

for line in open(sys.argv[1]):
    fname = line.split()[-1]
    rmsds = open(fname).readlines()
    if len(rmsds) == 0:
        print(line.rstrip())
    elif '_s5_' in line and len(rmsds) < 37:
        #this can legitimately happen if fewer than 9 poses are sampled per a run
        print(line.rstrip())

