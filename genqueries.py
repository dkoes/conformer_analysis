#!/usr/bin/env python3

import json, sys, os
from itertools import combinations

fname = sys.argv[1]

prefix,ext = os.path.splitext(fname)
q = json.load(open(fname))
features = []

for feat in  q['points']:
    if feat['enabled']:
        feat['radius'] = 1.0
        if 'size' in feat:
            del feat['size']
        if 'svector' in feat:
            del feat['svector']
        if 'vectr' in feat:
            del feat['vector']
        features.append(feat)

for i in range(3,len(features)+1):
    for j,combo in enumerate(combinations(features,i)):
        query = {"points": combo}
        json.dump(query, open(f'{prefix}_{i}_{j}.json','w'))

