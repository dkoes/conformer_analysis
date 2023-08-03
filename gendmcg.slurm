#!/bin/bash
#SBATCH --job gendmcg
#SBATCH --partition=any_gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:1

hostname
cd $SLURM_SUBMIT_DIR

rsync -av checkpoint_94.pt /tmp/
cmd=`sed -n "${SLURM_ARRAY_TASK_ID}p" gendmcgunmin`
eval $cmd

