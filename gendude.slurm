#!/bin/bash
#SBATCH --job gendude
#SBATCH --partition=dept_cpu
#SBATCH --nodes=1
#SBATCH -c 16 
##SBATCH --gres=gpu:1


cd $SLURM_SUBMIT_DIR

cmd=`sed -n "${SLURM_ARRAY_TASK_ID}p" gendude`
eval $cmd

