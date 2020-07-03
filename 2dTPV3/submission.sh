#!/bin/bash

# Job Name and Files (also --job-name)

#SBATCH -J slip4
#Output and error (also --output, --error):
#SBATCH -o ./job.%j.%x.out
#SBATCH -e ./job.%j.%x.err

#Initial working directory (also --chdir):
#SBATCH --workdir=/hppfs/work/pr63qo/di52lak2/ExaHyPE2020/ApplicationExamples/GPRDR/GPRDR_slip4/

#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=dli@geophysik.uni-muenchen.de

# Wall clock limit:
#SBATCH --time=00:30:00
#SBATCH --no-requeue

#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --account=pr63qo
#constraints are optional
#--constraint="scratch&work"
#SBATCH --partition=test

#Number of nodes and MPI tasks per node:
#max33 so far, else error
#SBATCH --nodes=14
#SBATCH --ntasks=83
##SBATCH --ntasks-per-node=6
## #SBATCH --cpus-per-task=8

#Needs specific MPI
#module switch mpi.intel mpi.intel/2019
#Run the program:

module load slurm_setup
module load mpi.intel/2019 tbb/2019 gcc/4.9

source /etc/profile.d/modules.sh
export MP_SINGLE_THREAD=no

srun --export=ALL ./ExaHyPE-GPRDR ./test2.exahype
