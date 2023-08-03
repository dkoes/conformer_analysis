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
Note that 77 structures failed to be parsed with RDKit.  We ignore these (newer version of RDKit are less problematic).

**COD (Crystallography Open Database)**
Original source is http://www.crystallography.net/cod/
Processed files (from https://chemrxiv.org/engage/chemrxiv/article-details/6356b98e2e0c63bd4d3a35f4) can be downloaded:
```
wget http://bits.csb.pitt.edu/files/cod.tgz
mkdir cod
cd cod
tar xvfz ../cod.tgz
```


**DUDE**
Downloaded from https://dude.docking.org/
```
wget https://dude.docking.org/db/subsets/all/all.tar.gz
tar xvfz all.tar.gz
mv all DUDE
```

**Crossdocking benchmark**
Downloaded from (original source  https://doi.org/10.1002/pro.3784)
```bash
wget https://bits.csb.pitt.edu/files/gnina1.0_paper/crossdocked_ds_data.tar.gz
tar xvfz crossdocked_ds_data.tar.gz
./make_wierbowski_smiles.py
curl https://raw.githubusercontent.com/gnina/models/master/data/Gnina1.0/ds_cd_input_pairs.txt | sed 's/carlos_cd/wierbowski_cd/g' | sed 's/PDB_Structures\///g' > ds_cd_input_pairs.txt 
```

## RDKit Conformer Generation

First we randomly sample a max of 250 conformers for each smile.  The commandlines can be generated with:
```bash
for f in platinum2017/*.smi refined-set/*/*.smi  wierbowski_cd/*/*.smi cod/*/*.smi
do 
 echo ./gen_rdkit_conformers.py --smi $f
done
```
If you are very patient you can remove the echo, or alternatively divide up the commands using your favorite batch processing technique (we use the `genrd.slurm` script to distribute across our SLURM cluster).

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
for f in platinum2017/*.smi refined-set/*/*.smi wierbowski_cd/*/*.smi cod/*/*.smi
do 
 echo ./dmcg.py --smi $f  --out ${f%.smi}_dmcg_250_unmin.sdf.gz --maxconfs 250 --eval-from /tmp/checkpoint_94.pt
done
```
This script also creates a `.txt` file with the RMSDs to the reference calculated by obrmsd.


## Energy Minimization

We next energy minimize the generated conformers using the UFF forcefield.  The resulting `_min.sdf.gz` file
is unordered.

```bash
for f in platinum2017/*unmin*sdf.gz refined-set/*/*unmin*.sdf.gz wierbowski_cd/*/*unmin*.sdf.gz cod*/*/*unmin.sdf.gz
do 
 echo ./minimize_sdf.py --sdf $f
done 
```
RMSDs to the reference conformation (if available) and all cross RMSDS (the RMSD between every pose in the file with every other pose - this is used for filtering) are also calculated with obrms and put into `.txt` and `.npy` files.  Note that a few cod conformers take a huge amount of time to calculate RMSDs for because of their symmetries.



## Bioactive Conformation Analysis.

See the notebook bioactive.ipynb.
We evaluate the effect the size of the conformational ensemble has on bioactive conformation identification as well as various methods for constructing ensembles (e.g., energy minimization and RMSD filtering).

## Pharmacophore Matching

### DUDE Conformer Generation

The `gen_dude_conformers.py` script will generate the following conformers:
 * Minimized but unsorted: 1,5,10,25,50,75,100,200
 * Minimized, sorted filters (0.5, 1.0, 2.0 RMSD) 1,5,10,25,50,75,100,200


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
Note we set each feature radius to 1.0 and remove any other constraints on the features.

For all identified interacting features, we then enumerate every possible
combination with at least three features.  Rather than worrying about identifying
a single best pharmacophore, we will simply test them all and take the best result.

```bash
for i in DUDE/*
do 
  pharmit pharma -receptor $i/receptor.pdb -in $i/crystal_ligand.mol2 -out $i/pharmit.json
  ./genqueries.py $i/pharmit.json
done
```
**Note:** Using this procedure only two interaction features were identified for the cah2 target, so the remaining two non-interacting features were added to ensure there are at least 3 to choose from.

### Pharmacophore Search

Setup parameters for building databases.

```bash
cd DUDE
for d in * ; do for i in 1 5 10 25 50 100 200; do echo "$d $i filtered_2.0"; echo "$d $i filtered_1.5"; echo "$d $i filtered_1.0"; echo "$d $i filtered_0.5"; echo "$d $i unranked"; done; done > params
mv params ../
cd ..
```

Count the number of conformers (sdsorter is available here: https://sourceforge.net/projects/sdsorter/)
```bash
for sdf in DUDE/*/*.sdf.gz
do
echo $sdf
sdsorter -printCnt $sdf > ${sdf%.sdf.gz}.cnt
sdsorter -reduceconfs 1 -printCnt $sdf >> ${sdf%.sdf.gz}.cnt
done
```

Build and search in parallel across a cluster.

`sbatch -a 1-3605 get_pharma_res.slurm`

### Analysis

See pharmacophore_analysis.ipynb


## Docking

First generate input conformations.  We will limit ourselves to a maximum ensemble 
size of 5 due to the cost of docking.  We generate these small ensembles:
 * 1 unminimized
 * 1 unbiased sampled minimized (same as above)
 * 1 minimal energy
 * 5 unbiased sampled energy minimized
 * 5 sorted by energy and filtered by RMSDs: 0, 0.5, 1.0, 2.0

for i in  refined-set/*/*unmin.sdf.gz wierbowski_cd/*/*unmin.sdf.gz 
do 
 echo ./setup_docking_confs.py $i; 
done > makedock

Remove water and other hetatoms from receptor structures of PDBbind.

```bash
for i in refined-set/*/*_protein.pdb
 do
 grep -v HETATM $i > ${i%_protein.pdb}_rec.pdb
done
```

Create a file containing all the docking commands:
```bash
./setup_docking_cmds.py > alldock
```

We run this using `dodock.slurm`.

To perform docking with smina, can modify the command lines as follows:
```bash
sed 's/gnina/smina/g' alldock | sed 's/_docked/_sdocked/g' > salldock
```


### Analysis

See docking_analysis.ipynb
