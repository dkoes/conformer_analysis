#!/bin/bash
#SBATCH --job dooode 
#SBATCH --partition=dept_cpu
#SBATCH --nodes=1

echo Running on `hostname`
echo workdir $SLURM_SUBMIT_DIR
echo ld_library_path $LD_LIBRARY_PATH

#scratch drive folder to work in
SCRDIR=/scr/${SLURM_JOB_ID}
#if the scratch drive doesn't exist (it shouldn't) make it.
if [[ ! -e $SCRDIR ]]; then
        mkdir $SCRDIR
fi
chmod +rX $SCRDIR
echo scratch drive ${SCRDIR}


cd $SLURM_SUBMIT_DIR


read -r dir i res <<< `sed -n ${SLURM_ARRAY_TASK_ID}p params`

echo "$d - $i - $res"

cp -r DUDE/$dir $SCRDIR/$dir
cd $SCRDIR/$dir


gunzip actives_final_rdkit_250_${res}.sdf.gz 
gunzip decoys_final_rdkit_250_${res}.sdf.gz 


echo "Building $i"
rm -rf db_${res}_${i}/ ; pharmit dbcreate  -in decoys_final_rdkit_250_${res}.sdf -in actives_final_rdkit_250_${res}.sdf -reduceconfs=${i} -dbdir db_${res}_${i}
echo "Searching $i"
for f in *.json; do echo -n "$f $i "; $SLURM_SUBMIT_DIR/getf1.py $f db_${res}_${i}; done > result_${res}_${i}.txt

cp result_${res}_${i}.txt $SLURM_SUBMIT_DIR/DUDE/$dir

sort -n -k4 result_${res}_${i}.txt  | tail -1

rm -rf db_* *.sdf *.gz

