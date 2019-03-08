#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J helium

#Output and error (also --output, --error):
#SBATCH -o ./%j.%x.out
#SBATCH -e ./%j.%x.err

#Initial working directory (also --chdir):
#SBATCH -D ./

#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=josef.winter@tum.de

# Wall clock limit:
#SBATCH --time=12:00:00
#SBATCH --no-requeue

#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn56re

#Partition
#SBATCH --partition=micro

#Number of nodes and MPI tasks per node:
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
 
#Important modules
module unload mkl mpi.intel intel
module load intel/19.0
module load gcc/9
module load mpi.intel/2019_gcc
module load cmake
module load slurm_setup
module load hdf5/1.8.20-intel-impi-threadsafe

export ALPACA_ENVIRONMENT=SNG
export COMPILER=GCC

#Run the program:
mpiexec -n $SLURM_NTASKS ./ALPACA Your/Input/File.xml 


