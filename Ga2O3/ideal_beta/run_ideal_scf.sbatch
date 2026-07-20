#!/bin/bash
#SBATCH --job-name=BetaSCF
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
mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np $SLURM_NTASKS /shared/quantum-espresso/bin/pw.x \
  -in beta_Ga2O3.scf_ideal.in > beta_Ga2O3.scf_ideal.out

date
