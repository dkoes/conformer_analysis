#!/usr/bin/env python3
'''Given an sdf file, UFF minimize all the conformations in it.  Save the result
in the same order with an energy property on each conformation.
'''

import numpy as np
import os, argparse,glob, sys
from rdkit.Chem import AllChem as Chem
import gzip, re, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--sdf", type=str, required=True,help="sdf file")

args = parser.parse_args()
fname = args.sdf

m = re.search(r'(.*?(_ligand)?)_',fname)
refsdf = m.group(1)+'.sdf'
refsdf = refsdf.replace('_nosc','')

out = fname.replace('.sdf.gz','.rmsds.txt')
subprocess.check_output(f'obrms -f -m {refsdf} {fname} > {out}',shell=True)

#cross rmsds
out = fname.replace('sdf.gz','crms.npy')
text = subprocess.check_output(f'obrms -m -x {fname}',shell=True)     
A = np.array([line.split(',')[1:] for line in text.decode().split('\n') if line],dtype=float)   
np.save(out,A)

