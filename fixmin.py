#!/usr/bin/env python3


import numpy as np
import os, argparse,glob, sys
from rdkit.Chem import AllChem as Chem
import gzip, re, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--sdf", type=str, required=True,help="sdf file")

args = parser.parse_args()

fname = args.sdf
out = fname.replace('sdf.gz','crms.npy')
text = subprocess.check_output(f'obrms -m -x {fname}',shell=True)     
A = np.array([line.split(',')[1:] for line in text.decode().split('\n') if line],dtype=float)   
np.save(out,A)    
