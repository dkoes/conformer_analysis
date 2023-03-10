#!/usr/bin/env python3

import numpy as np
import os, argparse,glob, sys
from rdkit import Chem
import multiprocessing, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--pdbbind_dir", type=str, default="refined-set",help="Unpacked directory of PDBbind Refined")
args = parser.parse_args()

def makesmi(fname):
    try:
        base = os.path.splitext(fname)[0]

        if fname.endswith('mol2'):
            mol = Chem.MolFromMol2File(fname)
        else:
            mol = next(Chem.SDMolSupplier(fname))


        out = open(base+'.smi','wt')

        out.write(Chem.MolToSmiles(mol))
        out.write('\t'+base.split('/')[-1])
        out.close()
        
        Chem.RemoveStereochemistry(mol)
        out = open(base+'_nosc.smi','wt')
        out.write(Chem.MolToSmiles(mol))
        out.write('\t'+base.split('/')[-1])
        out.close()
        
    except:
        print("Problem with",fname)
        return fname
        
pool = multiprocessing.Pool();
failures = list(filter(lambda x: x, pool.map(makesmi,glob.glob(f'{args.pdbbind_dir}/*/*ligand.mol2'))));        
print(f'{len(failures)} errors processing mol2 files')
morefailures = []
for fname in failures:
    fname = os.path.splitext(fname)[0]+'.sdf'
    if makesmi(fname) == fname:
        print("Still problematic:",fname)
        morefailures.append(fname)
print(f'{len(morefailures)} unsolvable failures')
