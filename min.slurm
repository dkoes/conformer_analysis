#!/bin/bash
#SBATCH --job minrd
#SBATCH --partition=any_cpu
#SBATCH --nodes=1


cd $SLURM_SUBMIT_DIR

cmd=`sed -n "${SLURM_ARRAY_TASK_ID}p" minsdf`
eval $cmd

