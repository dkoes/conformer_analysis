#!/bin/bash
#SBATCH --job dock 
#SBATCH --partition=any_cpu
#SBATCH -c 2
#SBATCH --nodes=1

export PATH=/net/pulsar/home/koes/dkoes/git/smina/build/:$PATH

cd $SLURM_SUBMIT_DIR

cmd=`sed -n "${SLURM_ARRAY_TASK_ID}p" alldock`
eval $cmd

