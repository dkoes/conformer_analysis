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

## RDKit Conformer Generation

First we randomly sample a max of 200 conformers for each smile.  The commandlines can be generated with:
```bash
for f in platinum2017/*.smi refined-set/*/*.smi 
do 
 echo ./gen_rdkit_conformers.py --smi $f
done
```
If you are very patient you can remove the echo, or alternative divide up the commands using your favorite batch processing technique (we use the genrd.slurm script to distribute across our SLURM cluster).

Request max 200, unminimized

Minimize the 200 

Create sorted subsets  1,10,25,50,100 with and without RMSD filtering: 0.25, 0.5, 0.75, 1.0

Create subsets of both minimized and unminimized: 1,10,25,50,100

## DMCG Conformer Generation

Checkout code and install

Use dmcg script to generate 200

Minimize and sort and make subsets


## Bioactive 

