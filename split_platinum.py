#!/usr/bin/env python3

import numpy as np
import os, argparse
from rdkit import Chem

parser = argparse.ArgumentParser()
parser.add_argument("--sdf", type=str, required=True,help="Input Platinum SDF")
parser.add_argument("--outdir",type=str, default="platinum2017",help="Output directory")
args = parser.parse_args()

os.makedirs(args.outdir,exist_ok=True)

for mol in Chem.SDMolSupplier(args.sdf):
    name = mol.GetProp('_Name')
    #the following is needed for chiral centers to be included in smiles output
    Chem.FindMolChiralCenters(mol,includeUnassigned=True)
    _,pdb,ch = name.split('_')
   #write sdf
    writer = Chem.SDWriter(f'{args.outdir}/{pdb}.sdf')
    writer.write(mol)
    writer.close()
    #write smi
    smi = Chem.MolToSmiles(mol)
    out = open(f'{args.outdir}/{pdb}.smi','wt')
    out.write(f'{smi} {pdb}\n')
    out.close()

    Chem.RemoveStereochemistry(mol)
    out = open(f'{args.outdir}/{pdb}_nosc.smi','wt')
    out.write(Chem.MolToSmiles(mol))
    out.write(f' {pdb}\n')
    out.close()    
