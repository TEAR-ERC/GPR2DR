#!/bin/bash


for i in  30 31 32
do
       echo $i
	head -310 GPRmaterials.f90.bk > GPRmaterials.f90
	echo "          RuptureCoeff(7) = "$i"         ! alpha2 " >> GPRmaterials.f90
	tail -294 GPRmaterials.f90.bk >> GPRmaterials.f90 

       if [ ! -d output$i ];then
        mkdir output$i
       fi

	head -155 test2.exahype > output${i}/test3.exahype
        echo '	"output": "./output'$i'/conserved",' >> output${i}/test3.exahype
        tail -16 test2.exahype >> output${i}/test3.exahype
      
 	../../../Toolkit/toolkit.sh output${i}/test3.exahype
        make -j32
        mv ExaHyPE-GPRDR output${i}/

    cat <<EOF > output${i}/submission2.sh
#!/bin/bash

# Job Name and Files (also --job-name)

#SBATCH -J Exahype 
#Output and error (also --output, --error):
#SBATCH -o ./job.%j.%x.out
#SBATCH -e ./job.%j.%x.err

#Initial working directory (also --chdir):
#SBATCH --workdir=/hppfs/work/pr63qo/di52lak2/ExaHyPE2020/ApplicationExamples/GPRDR/GPRDR/

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
module load mpi.intel/2019 tbb
source /etc/profile.d/modules.sh
export MP_SINGLE_THREAD=no
EOF

   echo "srun --export=ALL ./output"$i"/ExaHyPE-GPRDR ./output"$i"/test3.exahype" >> output$i/submission2.sh

        sbatch output$i/submission2.sh
done
