#!/usr/bin/env python3
'''Given a smile, generate minimized RDKIT conformers using ETKDGv3.
 Possibly rank by energy and filter by rmsd.
'''

import numpy as np
import os, argparse,glob, sys
from rdkit.Chem import AllChem as Chem
import gzip, subprocess, multiprocessing

parser = argparse.ArgumentParser()
parser.add_argument("--smi", type=str, required=True,help="smile file")
parser.add_argument("--maxconfs",type=int, default=250,help="max number of conformers to generate")
parser.add_argument("--seed",type=int,default=6262009,help="random seed")
parser.add_argument('--j',type=int,default=8,help='number of threads')
parser.add_argument('--filter',type=float,default=0,help='RMSD to filter by; if zero no sorting or ranking')
args = parser.parse_args()

base,_ = os.path.splitext(args.smi)

if args.filter == 0:
    ename = f'{base}_rdkit_{args.maxconfs}_unranked.sdf.gz'
else:
    ename = f'{base}_rdkit_{args.maxconfs}_filtered_{args.filter}.sdf.gz'
    
outetkdg = gzip.open(ename,'wt')
sdwriteretkdg = Chem.SDWriter(outetkdg)

def getRMS(mol, c1,c2):
    rms = Chem.GetBestRMS(mol,mol,c1,c2)
    return rms    
    
def create_confs(line):    
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
                
    #etkdg
    params = Chem.ETKDGv3()
    params.randomSeed = args.seed
    cids = Chem.EmbedMultipleConfs(mol, args.maxconfs, params)
    cenergy = []
    for conf in cids: 
        Chem.UFFOptimizeMolecule(mol,confId=conf)
        cenergy.append(Chem.UFFGetMoleculeForceField(mol,confId=conf).CalcEnergy())
        
    mol = Chem.RemoveHs(mol)
    
    mol.SetProp("_Name",name);
    
    sortedcids = sorted(cids,key = lambda cid: cenergy[cid])
    if len(sortedcids) > 0:
        mine = cenergy[sortedcids[0]]
    else:
        mine = 0
    if(args.filter > 0):
        selected = set()
        energies = []
        newmol = Chem.Mol(mol)
        newmol.RemoveAllConformers()
        for conf in sortedcids:
            #check rmsd
            passed = True
            for seenconf in selected:
                rms = getRMS(mol,seenconf,conf) 
                if rms < args.filter:
                    break
            else:
                selected.add(conf)
                energies.append(cenergy[conf])
                newmol.AddConformer(mol.GetConformer(conf),assignId=True)
        return newmol,name,energies
    else:
        return mol,name,cenergy
    
    



pool = multiprocessing.Pool(args.j)

for mol,name,energies in pool.imap_unordered(create_confs, open(args.smi,'rt')):
    #I have no idea what the name property doesn't stick through pickling..
    mol.SetProp('_Name',name)
    for cid in range(mol.GetNumConformers()):
        mol.SetDoubleProp('energy',energies[cid])
        sdwriteretkdg.write(mol, confId=cid)    


sdwriteretkdg.close()
outetkdg.close()

    
