#!/usr/bin/env python3

import numpy as np
import os, argparse, gzip
from rdkit import Chem

parser = argparse.ArgumentParser()
parser.add_argument("--sdf", type=str, required=True,help="Input  SDF")
parser.add_argument("--outdir",type=str, default="refined-set",help="Output directory")
args = parser.parse_args()

os.makedirs(args.outdir,exist_ok=True)

if args.sdf.endswith('.gz'):
    infile = gzip.open(args.sdf)
else:
    infile = open(args.sdf)

oldname = None
writer = None
outgz = None
for mol in Chem.ForwardSDMolSupplier(infile):
    name = mol.GetProp('_Name')
    if name != oldname:        
        oldname = name
        #create sdf
        if writer:
            writer.close()
            outgz.close()
        outgz = gzip.open(f'{args.outdir}/{name}/{name}_ligand_auto3d_250_min.sdf.gz','wt')
        writer = Chem.SDWriter(outgz)
    writer.write(mol)
    
 
