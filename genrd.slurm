#!/bin/bash
#SBATCH --job genrd
#SBATCH --partition=any_cpu
#SBATCH --nodes=1
##SBATCH --gres=gpu:1


cd $SLURM_SUBMIT_DIR

cmd=`sed -n "${SLURM_ARRAY_TASK_ID}p" genrdunmin`
eval $cmd

