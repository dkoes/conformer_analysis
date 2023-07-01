#!/usr/bin/env python3

'''Given a full ensemble of unminimized conformers, generate
smaller ensembles:
 1 unminimized
 1 unbiased sampled minimized (same as above)
 1 minimal energy
 5 identical copies of the minal energy
 5 unbiased sampled energy minimized
 5 sorted by energy and filtered by RMSDs: 0.5, 1.0, 2.0
 '''
import argparse
import glob, sys, os, re, gzip, io
from rdkit.Chem import AllChem as Chem
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("unmin", type=str, help="Unminimized filename; naming conventions are assumed for finding other files")
args = parser.parse_args()

m = re.search(r'(.*)_unmin.sdf.gz',args.unmin)
if not m:
    sys.stderr.write(f'{args.unmin} is not in the expected form (*_unmin.sf.gz)\n')
    sys.exit(-1)
base = m.group(1)

minfile = f'{base}_min.sdf.gz'
rmsfile = f'{base}_min.crms.npy'
ename = f'{base}_min.energy.npy'

energies = np.load(ename)
R = np.load(rmsfile)


def write_mols(fname, mols):
    '''Write mols to a gzip sdf'''
    out = gzip.open(fname,'wt')
    sdwriter = Chem.SDWriter(out)
    for mol in mols:
        sdwriter.write(mol)
    sdwriter.close()
    out.close()

#A single unmin
supp = Chem.ForwardSDMolSupplier(gzip.open(args.unmin))
mol = next(supp)
write_mols(f'{base}_2dock_single_unranked_unmin.sdf.gz',[mol])

#read all minimized
supp = Chem.ForwardSDMolSupplier(gzip.open(minfile))
minmols = [mol for mol in supp]

#The above minimized
write_mols(f'{base}_2dock_single_unranked_min.sdf.gz', minmols[:1])

#The lowest energy
e = np.load(ename)
minmol = minmols[np.argmin(e)] 
write_mols(f'{base}_2dock_single_eranked_min.sdf.gz', [minmol])

#The lowest energy duplicated 5 times
write_mols(f'{base}_2dock_5_same_min.sdf.gz', [minmol]*5)

#5 random min
write_mols(f'{base}_2dock_5_unranked_min.sdf.gz', minmols[:5])


def filter_poses(mols, energy, threshold, R, N):
    '''Take the N first lowest energies that are more than threshold apart'''
    sindex = np.argsort(energy)
    selected = [] #indices

    for i in sindex:
        for s in selected:
            if R[i,s] < threshold:
                break
        else: # did not break from loop
            selected.append(i)
        if len(selected) == N:
            break
    return np.array(mols)[selected]
    
    
for t in [0,0.5,1.0,2.0]:
    filtered = filter_poses(minmols, e, t, R, 5)
    write_mols(f'{base}_2dock_5_filtered_{t}_min.sdf.gz', filtered)
