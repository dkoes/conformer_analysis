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

base,_ = os.path.splitext(args.sdf)
if base.endswith('.sdf'):
    base,_ = os.path.splitext(base)

if '_unmin' in base:
    fname = base.replace('_unmin','_min')+'.sdf.gz'
else:
    fname = base+'_min,sdf.gz'
    
    
out=gzip.open(fname,'wt')
sdwriter = Chem.SDWriter(out)

if args.sdf.endswith('.gz'):
    infile = gzip.open(args.sdf)
else:
    infile = open(args.sdf)
supp = Chem.ForwardSDMolSupplier(infile)

energies = []
for mol in supp:
    try:
        mol = Chem.AddHs(mol,addCoords=True)
        Chem.UFFOptimizeMolecule(mol)
        energy = Chem.UFFGetMoleculeForceField(mol).CalcEnergy()
        mol.SetProp('energy',str(energy))
        mol = Chem.RemoveHs(mol)
        energies.append(energy)        
        sdwriter.write(mol)
    except Exception as e:
        print(args.sdf,e)

sdwriter.close()
out.close()

if 'wierbowski_cd' in fname:
    m = re.search(r'(.*?_LIG_aligned)',fname)
else:
    m = re.search(r'(.*?(_ligand)?)_',fname)

refsdf = m.group(1)+'.sdf'
refsdf = refsdf.replace('_nosc','')

out = fname.replace('.sdf.gz','.rmsds.txt')
subprocess.check_output(f'obrms -t 60 -f -m {refsdf} {fname} > {out}',shell=True)

#cross rmsds
out = fname.replace('sdf.gz','crms.npy')
text = subprocess.check_output(f'obrms -t 1 -m -x {fname}',shell=True)     
A = np.array([line.split(',')[1:] for line in text.decode().split('\n') if line],dtype=float)   
np.save(out,A)

#save energies
out = fname.replace('sdf.gz','energy.npy')
np.save(out, np.array(energies))
