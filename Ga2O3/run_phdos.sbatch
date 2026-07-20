#!/bin/bash
#SBATCH --job-name=ga2o3_phdos
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --output=%j_phdos.log

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
ulimit -s unlimited

source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl

mpirun -np "${SLURM_NTASKS}" /shared/quantum-espresso/bin/matdyn.x \
  -in Ga2O3.phdos.in > Ga2O3.phdos.out
python3 plot_ga2o3_phdos.py
