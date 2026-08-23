#!/bin/bash
#SBATCH --job-name=Fe2O3-BMeq-magcheck
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --array=0-1
#SBATCH --output=%A_%a.log

set -euo pipefail
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
STATES=(AFM FM)
STATE=${STATES[$SLURM_ARRAY_TASK_ID]}
cd "$SLURM_SUBMIT_DIR"
mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/pw.x \
  -npool 2 -ndiag 1 -in "$STATE.in" > "$STATE.out"
grep -q 'convergence has been achieved' "$STATE.out"
grep -q 'JOB DONE' "$STATE.out"
