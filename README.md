# Analysis of Conformer Generation Procedures for Structure Based Drug Discovery

It is assumed RDKit is installed (e.g. `pip install rdkit-pypi`).  We used version 2022.9.

## Raw Data

**Platinum 2017 dataset**, downloaded March 10, 2023
```
wget https://comp3d.univie.ac.at/fileadmin/user_upload/p_comp3d/datasets/platinum-dataset-2017-01-sdf.zip
unzip platinum-dataset-2017-01-sdf.zip
```
Separate into indvidual files and generate smiles:
```
./split_platinum.py --sdf platinum_dataset_2017_01.sdf
```

**PDBbind 2020**
You will need to register and download the refined tarball at [http://www.pdbbind.org.cn/download.php](http://www.pdbbind.org.cn/download.php)
```
tar xvfz PDBbind_v2020_refined.tar.gz 
./make_refined_smiles.py 
```
Note that 77 structures failed to be parsed with RDKit.  We ignore these.

**DUDE**
Downloaded from https://dude.docking.org/
```
wget https://dude.docking.org/db/subsets/all/all.tar.gz
tar xvfz all.tar.gz
mv all DUDE
```

## RDKit Conformer Generation

First we randomly sample a max of 250 conformers for each smile.  The commandlines can be generated with:
```bash
for f in platinum2017/*.smi refined-set/*/*.smi 
do 
 echo ./gen_rdkit_conformers.py --smi $f
done
```
If you are very patient you can remove the echo, or alternative divide up the commands using your favorite batch processing technique (we use the `genrd.slurm` script to distribute across our SLURM cluster).

This results in files named `4EB8_rdkit_250_dg_unmin.sdf.gz` and `4EB8_rdkit_250_etkdg_unmin.sdf.gz`.  
This script also uses the `obrms` utility from OpenBabel to calculate the RMSD of each generated conformer and puts the result in `.txt` file.





## DMCG Conformer Generation

Checkout code and install (we are using commit 23e90f7231eac9a065093ee44b378cbf65181a7c).
```
git clone https://github.com/DirectMolecularConfGen/DMCG.git
cd DMCG
pip install .
```

You need to [download the Large Drugs model weights](https://drive.google.com/drive/folders/1piz0gy24bSt_e0ICjEcV5olFcCaJydwG).  This file is called `checkpoint_94.pt`.

DMCG has several dependencies (e.g. `torch_geometric`, `torch_sparse`, `torch_scatter` and [DGL](https://www.dgl.ai/pages/start.html)) that may need to be installed.  Make sure you can run the `dmcg.py` script on a test smile:
```
dmcg.py --smi test.smi --out test.sdf.gz --maxconfs 25 --eval-from checkpoint_94.pt
```

Setup generation commands the same as with rdkit:
```bash
for f in platinum2017/*.smi refined-set/*/*.smi 
do 
 echo ./dmcg.py --smi $f  --out ${f%.smi}_dmcg_250_unmin.sdf.gz --maxconfs 250 --eval-from /tmp/checkpoint_94.pt
done
```
This script also creates a `.txt` file with the RMSDs to the reference calculated by obrmsd.


## Energy Minimization

We next energy minimize the generated conformers using the UFF forcefield.  The resulting `_min.sdf.gz` file
is unordered.

```bash
for f in platinum2017/*unmin*sdf.gz refined-set/*/*unmin*.sdf.gz 
do 
 echo ./minimize_sdf.py --sdf $f
done 
```
RMSDs to the reference conformation and all cross RMSDS (the RMSD between every pose in the file with every other pose - this is used for filtering) are also calculated with obrms and put into `.txt` and `.npy` files.



## Bioactive Conformation Analysis.

See the notebook bioactive.ipynb.
We evaluate the effect the size of the conformational ensemble has on bioactive conformation identification as well as various methods for constructing ensembles (e.g., energy minimization and RMSD filtering).

## Pharmacophore Matching

### DUDE Conformer Generation
Minimized but unsorted: 1,5,10,25,50,75,100,200
Minimized, sorted filters (0.5, 1.0, 2.0 RMSD) 1,5,10,25,50,75,100,200
Energy window?

```bash
#label actives
for f in DUDE/*/actives_final.ism
do 
  sed -i 's/$/_active/' $f
done

for i in DUDE/*/*.ism
do 
  echo ./gen_dude_conformers.py --smi $i
done
```

### Pharmacophore Elucidation

Use Pharmit to identify the interacting features on the crystal ligand.
This uses SMARTS expressions to identify features and distance and count based
heuristics to enable them if they are interacting appropriately with the protein.

For all identified interacting features, when then enumerate every possible
combination with at least three features.  Rather than worrying about identifying
a single best pharmacophore, we will simply test them all and take the best result.

```bash
for i in DUDE/*
do 
  pharmit pharma -receptor $i/receptor.pdb -in $i/crystal_ligand.mol2 -out $i/pharmit.json
  ./genqueries.py $i/pharmit.json
done
```

### Pharmacophore Search

Setup parameters for building databases.

```bash
cd DUDE
for d in * ; do for i in 1 5 10 25 50 100 200; do echo "$d $i filtered_2.0"; echo "$d $i filtered_1.5"; echo "$d $i filtered_1.0"; echo "$d $i filtered_0.5"; echo "$d $i unranked"; done; done > params
mv params ../
cd ..
```

Build and search in parallel across a cluster.

`sbatch -a 1-3605 get_pharma_res.slurm`


