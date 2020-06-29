#!/bin/bash
#SBATCH --output=./slurm-%j.out
#SBATCH --error=./slurm-%j.err
#SBATCH --workdir=./
#SBATCH --job-name=slip2
#SBATCH --nodes=6
#SBATCH --ntasks=83
#SBATCH --mail-type=all
#SBATCH --time=24:00:00

#export OMP_NUM_THREADS=24
#export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS
export MP_SINGLE_THREAD=no
 
# run compute job
srun --export=ALL ./ExaHyPE-GPRDR test2.exahype


