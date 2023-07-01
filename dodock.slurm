#!/bin/bash
#SBATCH --job dock 
#SBATCH --partition=any_cpu
#SBATCH -c 1
#SBATCH --nodes=1


cd $SLURM_SUBMIT_DIR

cmd=`sed -n "${SLURM_ARRAY_TASK_ID}p" alldock`
eval $cmd

