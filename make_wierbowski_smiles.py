#!/usr/bin/env python3

import numpy as np
import os, argparse,glob, sys
from rdkit import Chem
import multiprocessing, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--dir", type=str, default="wierbowski_cd/",help="Unpacked directory of (downsampled) Wierbowski dataset")
args = parser.parse_args()

def makesmi(fname):
    try:
        base = os.path.splitext(fname)[0]
        mol = next(Chem.SDMolSupplier(fname))

        out = open(base+'.smi','wt')
        out.write(Chem.MolToSmiles(mol))
        out.write('\t'+base.split('/')[-1])
        out.close()
        
    except:
        print("Problem with",fname)
        return fname
        
pool = multiprocessing.Pool();
failures = list(filter(lambda x: x, pool.map(makesmi,glob.glob(f'{args.dir}/*/*LIG_aligned.sdf'))));        
print(f'{len(failures)} errors processing sdf files')

