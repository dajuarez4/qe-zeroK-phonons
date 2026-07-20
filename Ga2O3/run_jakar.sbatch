#!/bin/bash
#SBATCH --job-name=Phonns
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --output=%j.log

pwd
hostname
date
ulimit -s unlimited

# Load Modules
source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl

# Execute Program
#mpirun -np $SLURM_NTASKS /shared/quantum-espresso/bin/ph.x -in Ga2O3.ph.in > Ga2O3.ph.out
mpirun --mca pml ob1 --mca btl ^openib,uct -np $SLURM_NTASKS /shared/quantum-espresso/bin/matdyn.x -in Ga2O3.matdyn.in > Ga2O3.matdyn.out
