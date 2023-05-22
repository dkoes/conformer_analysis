#!/usr/bin/env python3
'''Given a smile, generate minimized RDKIT conformers using ETKDGv3.
 Rank by energy and then filter by different rmsd thresholds.
'''

import numpy as np
import os, argparse,glob, sys
from rdkit.Chem import AllChem as Chem
import gzip, subprocess, multiprocessing
import traceback

parser = argparse.ArgumentParser()
parser.add_argument("--smi", type=str, required=True,help="smile file")
parser.add_argument("--maxconfs",type=int, default=250,help="max number of conformers to generate")
parser.add_argument("--seed",type=int,default=6262009,help="random seed")
parser.add_argument('--j',type=int,default=32,help='number of threads')
args = parser.parse_args()

base,_ = os.path.splitext(args.smi)

thresholds = [0, 0.5, 1.0, 1.5, 2.0]
sdwriters = []
outs = []
for t in thresholds:  
  if t == 0:
      ename = f'{base}_rdkit_{args.maxconfs}_unranked.sdf.gz'
  else:
      ename = f'{base}_rdkit_{args.maxconfs}_filtered_{t}.sdf.gz'    
  outetkdg = gzip.open(ename,'wt')
  outs.append(outetkdg)
  sdwriteretkdg = Chem.SDWriter(outetkdg)
  sdwriters.append(sdwriteretkdg)

def getRMS(mol, c1,c2):
    rms = Chem.GetBestRMS(mol,mol,c1,c2)
    return rms    
    
def create_confs(line):   
 try:
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
    
    ret = {}
    for t in thresholds:
      if t > 0:
        selected = set()
        energies = []
        newmol = Chem.Mol(mol)
        newmol.RemoveAllConformers()
        for conf in sortedcids:
            #check rmsd
            passed = True
            for seenconf in selected:
                rms = getRMS(mol,seenconf,conf) 
                if rms < t:
                    break
            else:
                selected.add(conf)
                energies.append(cenergy[conf])
                newmol.AddConformer(mol.GetConformer(conf),assignId=True)
        ret[t] = (newmol,name,energies)
      else:
        ret[t] = (mol,name,cenergy)
    return ret
 except:
     print("Error",line)
     traceback.print_exc()
    



pool = multiprocessing.Pool(args.j)

nummols = {t:0 for t in thresholds}
numconfs = {t:0 for t in thresholds}

for res in pool.imap_unordered(create_confs, open(args.smi,'rt')):
    #I have no idea why the name property doesn't stick through pickling..
    if res == None:
        continue # error
    for t,writer in zip(thresholds,sdwriters):
      mol,name,energies = res[t]
      mol.SetProp('_Name',name)
      nummols[t] += 1
      for cid in range(mol.GetNumConformers()):
          mol.SetDoubleProp('energy',energies[cid])
          writer.write(mol, confId=cid)    
          numconfs[t] += 1


for writer in sdwriters:
  writer.close()
  
for t,out in zip(thresholds,outs):
  cntname = out.name.replace('.sdf.gz','.cnt')
  cout = open(cntname,'wt')
  cout.write(f'NumConfs: {numconfs[t]}\n')
  cout.write(f'NumMols: {nummols[t]}\n')
  cout.close()
  out.close()


    
