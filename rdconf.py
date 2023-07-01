#!/usr/bin/python3

import sys,string,argparse
from rdkit.Chem import AllChem as Chem
from optparse import OptionParser
import os, gzip

'''Given a smiles file, generate 3D conformers in output sdf.  
Energy minimizes and filters conformers to meet energy window and rms constraints.

Some time ago I compared this to alternative conformer generators and
it was quite competitive (especially after RDKit's UFF implementation 
added OOP terms).
'''

#convert smiles to sdf
def getRMS(mol, c1,c2):
    rms = Chem.GetBestRMS(mol,mol,c1,c2)
    return rms

parser = OptionParser(usage="Usage: %prog [options] <input>.smi <output>.sdf")
parser.add_option("--maxconfs", dest="maxconfs",action="store",
                  help="maximum number of conformers to generate per a molecule (default 20)", default="20", type="int", metavar="CNT")
parser.add_option("--sample_multiplier", dest="sample",action="store",
                  help="sample N*maxconfs conformers and choose the maxconformers with lowest energy (default 1)", default="1", type="float", metavar="N")
parser.add_option("--seed", dest="seed",action="store",
                  help="random seed (default 9162006)", default="9162006", type="int", metavar="s")
parser.add_option("--energy_window", dest="energy",action="store",
                  help="filter based on energy difference with lowest energy conformer", default="10", type="float", metavar="E")
parser.add_option("-v","--verbose", dest="verbose",action="store_true",default=False,
                  help="verbose output")             
parser.add_option("--nomin", dest="nomin",action="store_true",default=False,
                  help="don't perform energy minimization (bad idea)")                  


    
(options, args) = parser.parse_args()

if(len(args) < 2):
    parser.error("Need input and output")
    sys.exit(-1)
    
input = args[0]
output = args[1]
smifile = open(input)


split = os.path.splitext(output)
if split[1] == '.gz':
    outf=gzip.open(output,'wt+')
    output = split[0] #strip .gz
else:
    outf = open(output,'w+')
    
 
if os.path.splitext(output)[1] == '.pdb':
    sdwriter = Chem.PDBWriter(outf)
else:
    sdwriter = Chem.SDWriter(outf)
    
if sdwriter is None:
    print("Could not open ".output)
    sys.exit(-1)
    
for line in smifile:
    toks = line.split()
    smi = toks[0]
    name = ' '.join(toks[1:])    
    
    pieces = smi.split('.')
    if len(pieces) > 1:
        smi = max(pieces, key=len) #take largest component by length
        print("Taking largest component: %s\t%s" % (smi,name))
        
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        if options.verbose:
            print(smi)
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
            mol.SetProp("_Name",name);
            
            maxconfs = int(options.maxconfs)
            
            params = Chem.ETKDGv3()
            params.ignoreSmoothingFail = True

            for i in range(maxconfs):
                params.randomSeed = options.seed+i
                cids = Chem.EmbedMultipleConfs(mol, maxconfs, params)
                maxconfs -= len(cids)
                
                if options.verbose:
                    print(len(cids),"conformers found")
                cenergy = []            
                for conf in cids:
                    #not passing confID only minimizes the first conformer
                    if options.nomin:
                        cenergy.append(conf)
                    else:
                        converged = not Chem.UFFOptimizeMolecule(mol,confId=conf)
                        cenergy.append(Chem.UFFGetMoleculeForceField(mol,confId=conf).CalcEnergy())
                    if options.verbose:
                        print("Convergence of conformer",conf,converged,cenergy[-1])
                            
                for conf in cids:
                    sdwriter.write(mol,conf)
                if maxconfs == 0:
                    break
                    
        except (KeyboardInterrupt, SystemExit):
            raise                
        except Exception as e:
            print("Exception",e)
    else:
        print("ERROR:",smi)

sdwriter.close()
outf.close()
