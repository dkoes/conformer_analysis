#!/bin/bash
#SBATCH --job gendmcg
#SBATCH --partition=dept_gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH -x g005


cd $SLURM_SUBMIT_DIR

rsync -av checkpoint_94.pt /tmp/
cmd=`sed -n "${SLURM_ARRAY_TASK_ID}p" gendmcgunmin`
eval $cmd

