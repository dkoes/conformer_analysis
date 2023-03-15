#!/bin/bash
#SBATCH --job gendude
#SBATCH --partition=any_cpu
#SBATCH --nodes=1
#SBATCH -c 8
##SBATCH --gres=gpu:1


cd $SLURM_SUBMIT_DIR

cmd=`sed -n "${SLURM_ARRAY_TASK_ID}p" gendude`
eval $cmd

