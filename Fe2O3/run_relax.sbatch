#!/bin/bash
#SBATCH --job-name=Fe2O3Rel
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --output=%j.log

pwd
hostname
date
ulimit -s unlimited

source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np $SLURM_NTASKS /shared/quantum-espresso/bin/pw.x \
  -in alpha_Fe2O3.relax.in > alpha_Fe2O3.relax.out

if grep -q "JOB DONE" alpha_Fe2O3.relax.out; then
  python3 update_scf_from_relax.py
else
  echo "Relaxation did not finish cleanly; SCF input was not changed."
  exit 1
fi

date

