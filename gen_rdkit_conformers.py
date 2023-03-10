#!/usr/bin/env python3
'''Given a smile, generate unminimized RDKIT conformers using default distance geometry and ETKDG.
Files will be suffixed with _rdkit_200_dg_unmin.sdf.gz and _rdkit_200_etkdg_unmin.sdf.gz
'''

import numpy as np
import os, argparse,glob, sys
from rdkit.Chem import AllChem as Chem
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--smi", type=str, required=True,help="smile file")
parser.add_argument("--maxconfs",type=int, default=200,help="max number of conformers to generate")
parser.add_argument("--seed",type=int,default=6262009,help="random seed")
args = parser.parse_args()

base,_ = os.path.splitext(args.smi)

outdg=gzip.open(f'{base}_rdkit_{args.maxconfs}_dg_unmin.sdf.gz','wt+')
sdwriterdg = Chem.SDWriter(outdg)
outetkdg = gzip.open(f'{base}_rdkit_{args.maxconfs}_etkdg_unmin.sdf.gz','wt+')
sdwriteretkdg = Chem.SDWriter(outetkdg)

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
    mol2 = Chem.RemoveHs(mol2)
    for cid in range(mol2.GetNumConformers()):
        sdwriterdg.write(mol2, confId=cid)    
    
    #etkdg
    params = Chem.ETKDGv3()
    params.randomSeed = args.seed
    Chem.EmbedMultipleConfs(mol, args.maxconfs, params)
    mol = Chem.RemoveHs(mol)
    for cid in range(mol2.GetNumConformers()):
        sdwriteretkdg.write(mol2, confId=cid)    


sdwriterdg.close()
sdwriteretkdg.close()
