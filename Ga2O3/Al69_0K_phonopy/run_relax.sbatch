#!/bin/bash
#SBATCH --job-name=Al69-relax
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --output=%j.log

set -euo pipefail

pwd
hostname
date
ulimit -s unlimited

set +u
source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
set -u
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$SLURM_SUBMIT_DIR"
mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/pw.x \
  -npool 2 -ndiag 1 \
  -in 01_vc-relax.in > 01_vc-relax.out

grep -q 'JOB DONE' 01_vc-relax.out
date
