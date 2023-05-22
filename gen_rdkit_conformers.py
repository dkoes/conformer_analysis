#!/usr/bin/env python3
'''Given a smile, generate unminimized RDKIT conformers using default distance geometry and ETKDG.
Files will be suffixed with _rdkit_250_dg_unmin.sdf.gz and _rdkit_250_etkdg_unmin.sdf.gz
'''

import numpy as np
import os, argparse,glob, sys
from rdkit.Chem import AllChem as Chem
import gzip, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--smi", type=str, required=True,help="smile file")
parser.add_argument("--maxconfs",type=int, default=250,help="max number of conformers to generate")
parser.add_argument("--seed",type=int,default=6262009,help="random seed")
args = parser.parse_args()

base,_ = os.path.splitext(args.smi)

dgname = f'{base}_rdkit_{args.maxconfs}_dg_unmin.sdf.gz'
outdg=gzip.open(dgname,'wt')
sdwriterdg = Chem.SDWriter(outdg)

ename = f'{base}_rdkit_{args.maxconfs}_etkdg_unmin.sdf.gz'
outetkdg = gzip.open(ename,'wt')
sdwriteretkdg = Chem.SDWriter(outetkdg)

dgcnt = 0
ecnt = 0

for line in open(args.smi,'rt'):
    toks = line.split()
    smi = toks[0]
    name = ' '.join(toks[1:])    
    
    pieces = smi.split('.')
    if len(pieces) > 1:
        smi = max(pieces, key=len) #take largest component by length
        print("Taking largest component: %s\t%s" % (smi,name))
    
    mol = Chem.MolFromSmiles(smi)

    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    mol.SetProp("_Name",name);
                
    #distance geometry
    mol2 = Chem.Mol(mol)
    Chem.EmbedMultipleConfs(mol2, args.maxconfs,randomSeed=args.seed)
    if mol2.GetNumConformers() == 0:  #failed, try w/o sc
        Chem.RemoveStereochemistry(mol2)
        Chem.EmbedMultipleConfs(mol2, args.maxconfs,randomSeed=args.seed) # there is exactly one platinum molecule where this triggers    
    mol2 = Chem.RemoveHs(mol2)
    for cid in range(mol2.GetNumConformers()):
        sdwriterdg.write(mol2, confId=cid)    
        dgcnt += 1
    
    #etkdg
    params = Chem.ETKDGv3()
    params.randomSeed = args.seed
    Chem.EmbedMultipleConfs(mol, args.maxconfs, params)
    if mol.GetNumConformers() == 0:  #failed, try w/o sc
        Chem.RemoveStereochemistry(mol)
        Chem.EmbedMultipleConfs(mol, args.maxconfs,randomSeed=args.seed) # there is exactly one platinum molecule where this triggers        
    mol = Chem.RemoveHs(mol)
    for cid in range(mol.GetNumConformers()):
        sdwriteretkdg.write(mol, confId=cid)    
        ecnt += 1


sdwriterdg.close()
sdwriteretkdg.close()
outdg.close()
outetkdg.close()

if dgcnt == 0:
    os.remove(dgname)
if ecnt == 0:
    os.remove(ename)
    
#compute rmsds 
refsdf = base.replace('_nosc','')+'.sdf'
outdg = dgname.replace('.sdf.gz','.rmsds.txt')
text = subprocess.check_output(f'obrms -f -m {refsdf} {dgname} > {outdg}',shell=True)

oute = ename.replace('.sdf.gz','.rmsds.txt')
text = subprocess.check_output(f'obrms -f -m {refsdf} {ename} > {oute}',shell=True)
